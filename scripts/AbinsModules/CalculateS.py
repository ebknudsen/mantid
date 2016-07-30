import numpy as np

from IOmodule import IOmodule
from AbinsData import AbinsData
from QData import QData
from CalculateQ import CalculateQ
from InstrumentProducer import InstrumentProducer
from CalculatePowder import CalculatePowder
from CrystalData import CrystalData
from PowderData import  PowderData
from SData import SData
from Instruments import Instrument

import Constants

class CalculateS(IOmodule):
    """
    Class for calculating S(Q, omega)
    """

    def __init__(self, filename=None, temperature=None, sample_form=None, abins_data=None, instrument_name=None):
        """
        @param filename: name of input DFT file (CASTEP: foo.phonon)
        @param temperature: temperature in K for which calculation of S should be done
        @param sample_form: form in which experimental sample is: Powder or SingleCrystal (str)
        @param abins_data: object of type AbinsData with data from phonon file
        @param instrument_name: name of instrument (str)

        """
        super(CalculateS, self).__init__(input_filename=filename, group_name=Constants.S_data_group)

        if not (isinstance(temperature, float) or isinstance(temperature, int)):
            raise ValueError("Invalid value of the temperature. Number was expected.")
        if temperature < 0:
            raise ValueError("Temperature cannot be negative.")
        self._temperature = temperature

        if  sample_form in Constants.all_sample_forms:
            self._sample_form = sample_form
        else:
            raise  ValueError("Invalid sample form %s"%sample_form)

        if isinstance(abins_data, AbinsData):
            self._abins_data = abins_data
        else:
            raise ValueError("Object of type AbinsData was expected.")

        if instrument_name in Constants.all_instruments:
            self._instrument_name = instrument_name
        else:
            raise ValueError("Unknown instrument %s" % instrument_name)


    def _calculate_S(self):

        # Produce instrument object
        _instrument_producer = InstrumentProducer()
        _instrument = _instrument_producer.produceInstrument(name=self._instrument_name)

        # Calculate Q
        _q_calculator = CalculateQ(filename=self._input_filename, instrument=_instrument, sample_form=self._sample_form)
        if self._instrument_name != "None":

            _q_calculator.collectFrequencies(k_points_data=self._abins_data.getKpointsData())
        _q_vectors = _q_calculator.getQvectors()

        # Powder case: calculate MSD and DW
        if self._sample_form == "Powder":

            _powder_calculator = CalculatePowder(filename=self._input_filename, abins_data=self._abins_data, temperature=self._temperature)
            _powder_data = _powder_calculator.getPowder()
            _s_data = self._calculate_s_powder(q_data=_q_vectors, powder_data=_powder_data, instrument=_instrument)

        # Crystal case: calculate DW
        else:
            raise ValueError("SingleCrystal case not implemented yet.")

        return _s_data


    def _calculate_s_powder(self, q_data=None, powder_data=None, instrument=None):

        """
        Calculates S for the powder case.

        @param q_data:  data with Q vectors (or squared q vectors); it is an object of type QData
        @param powder_data: object of type PowderData with mean square displacements and Debye-Waller factors for the case of powder
        @param instrument: object of type instrument; instance of the instrument for which the whole simulation is performed
        @return: object of type SData with dynamical structure factor for the powder case
        """
        if not isinstance(q_data, QData):
            raise ValueError("Input parameter 'q_data'  should be of type QData." )

        if not isinstance(powder_data, PowderData):
            raise ValueError("Input parameter 'powder_data' should be of type PowderData.")

        if not isinstance(instrument, Instrument):
            raise ValueError("Input parameter 'instrument' should be of type Instrument.")

        _items = []

        _powder_atom_data = powder_data.extract()
        _num_atoms = _powder_atom_data["msd"].shape[0]
        _msd = _powder_atom_data["msd"]
        _dw = _powder_atom_data["dw"]

        _extracted_q_data = q_data.extract()

        _abins_data_extracted = self._abins_data.extract()
        _atom_data = _abins_data_extracted["atoms_data"]
        _k_points_data = _abins_data_extracted["k_points_data"]
        _num_k = _k_points_data["k_vectors"].shape[0]
        _num_freq = _k_points_data["frequencies"].shape[1]

        _factorials = [ np.math.factorial(n) for n in range(Constants.overtones_num)]

        _value = np.zeros( (Constants._pkt_per_peak * _num_freq, Constants.overtones_num + 1), dtype=Constants.float_type)
        _frequencies = np.zeros( (Constants._pkt_per_peak * _num_freq), dtype=Constants.float_type)



        _value_dft = np.zeros( _num_freq, dtype=Constants.float_type) # DFT discrete peaks
        _s_sum = Constants.overtones_num

        for atom in range(_num_atoms):

            _value.fill(0.0)

            for overtone in range(Constants.overtones_num):
                for k in range(_num_k):

                    # noinspection PyTypeChecker
                    # only discrete values
                    # comment:  (_extracted_q_data[k] * _dw[atom])^2 are equivalent to DW from  "Vibrational spectroscopy with neutrons...."
                    _value_dft = np.power(_extracted_q_data[k] * _msd[atom], 2 * overtone) * np.power(np.exp(- _extracted_q_data[k] * _dw[atom]), 2)

                    #  convolve value with instrumental resolution; resulting spectrum has broadened peaks with Gaussian-like shape
                    _value[:, overtone] += instrument.convolve_with_resolution_function(frequencies=_k_points_data["frequencies"][k],
                                                                                        s_dft=_value_dft,
                                                                                        points_per_peak=Constants._pkt_per_peak)



                _value[:, overtone] /=  _factorials[overtone] * 4 * np.pi # overtone contribution
                _value[:, _s_sum] += _value[:, overtone]  # total contribution

            _items.append({"sort":  _atom_data[atom]["sort"],
                           "symbol": _atom_data[atom]["symbol"],
                           "value":  np.copy(_value)})


        for k in range(_num_k):
            _frequencies += instrument.produce_abscissa(frequencies=_k_points_data["frequencies"][k], points_per_peak=Constants._pkt_per_peak) * _k_points_data["weights"][k]

        # Variable items is a list. Each entry in the list corresponds to one atom. Each  entry has a form of dictionary
        # with the following entries:
        #
        #     sort:    defines symmetry equivalent atoms (not functional at the moment)
        #     symbol:  symbol of atom (hydrogen -> H)
        #     value:   value  of dynamical structure factor for the given atom in powder scenario; first dimension
        #              stores frequencies and the second stores dynamical structure factor S
        #

        _s_data = SData(temperature=self._temperature, sample_form=self._sample_form)
        _s_data.set(items=dict(atoms=_items, frequencies=_frequencies))

        return _s_data


    def _calculate_s_crystal(self, crystal_data=None):

        if not isinstance(crystal_data, CrystalData):
            raise ValueError("Input parameter should be of type CrystalData.")
        # TODO: implement calculation of S for the single crystal scenario






    def getS(self):
        """
        Calculates dynamical structure factor S.
        @return: object of type SData.
        """

        data = self._calculate_S()

        self.addAttribute("temperature", self._temperature)
        self.addAttribute("sample_form", self._sample_form)
        self.addAttribute("filename", self._input_filename)
        extracted_data = data.extract()
        self.addStructuredDataset("atoms_data", extracted_data["atoms"])
        self.addNumpyDataset("convoluted_frequencies", extracted_data["frequencies"])

        self.save()

        return data


    def loadData(self):
        """
        Loads S from an hdf file.
        @return: object of type SData.
        """
        _data = self.load(list_of_structured_datasets=["atoms_data"],
                          list_of_numpy_datasets=["convoluted_frequencies"],
                          list_of_attributes=["temperature", "sample_form", "filename"])
        _s_data = SData(temperature=_data["attributes"]["temperature"],
                        sample_form=_data["attributes"]["sample_form"])
        _s_data.set(items=dict(atoms=_data["structured_datasets"]["atoms_data"],
                               frequencies=_data["datasets"]["convoluted_frequencies"]))

        return _s_data
