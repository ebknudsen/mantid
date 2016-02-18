#pylint: disable=invalid-name,relative-import,W0611,R0921,R0902,R0904,R0921,C0302
################################################################################
#
# MainWindow application for reducing HFIR 4-circle
#
################################################################################
import os
import math
import csv
import time
import random

from PyQt4 import QtCore, QtGui
from mantidqtpython import MantidQt

import reduce4circleControl as r4c
import guiutility as gutil
import fourcircle_utility as fcutil

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

# import line for the UI python class
from ui_MainWindow import Ui_MainWindow


class MainWindow(QtGui.QMainWindow):
    """ Class of Main Window (top)
    """
    def __init__(self, parent=None):
        """ Initialization and set up
        """
        # Base class
        QtGui.QMainWindow.__init__(self,parent)

        # UI Window (from Qt Designer)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # Make UI scrollable
        self._scrollbars = MantidQt.API.WidgetScrollbarDecorator(self)
        self._scrollbars.setEnabled(True) # Must follow after setupUi(self)!

        # Mantid configuration
        self._instrument = str(self.ui.comboBox_instrument.currentText())
        # config = ConfigService.Instance()
        # self._instrument = config["default.instrument"]

        # Event handling definitions
        # Top
        self.connect(self.ui.pushButton_setExp, QtCore.SIGNAL('clicked()'),
                     self.do_set_experiment)

        # Tab 'Data Access'
        self.connect(self.ui.pushButton_applySetup, QtCore.SIGNAL('clicked()'),
                     self.do_apply_setup)
        self.connect(self.ui.pushButton_browseLocalDataDir, QtCore.SIGNAL('clicked()'),
                     self.do_browse_local_spice_data)
        self.connect(self.ui.pushButton_testURLs, QtCore.SIGNAL('clicked()'),
                     self.do_test_url)
        self.connect(self.ui.pushButton_ListScans, QtCore.SIGNAL('clicked()'),
                     self.do_list_scans)
        self.connect(self.ui.pushButton_downloadExpData, QtCore.SIGNAL('clicked()'),
                     self.do_download_spice_data)
        self.connect(self.ui.comboBox_mode, QtCore.SIGNAL('currentIndexChanged(int)'),
                     self.change_data_access_mode)

        # Tab 'View Raw Data'
        self.connect(self.ui.pushButton_setScanInfo, QtCore.SIGNAL('clicked()'),
                     self.do_load_scan_info)
        self.connect(self.ui.pushButton_plotRawPt, QtCore.SIGNAL('clicked()'),
                     self.do_plot_pt_raw)
        self.connect(self.ui.pushButton_prevPtNumber, QtCore.SIGNAL('clicked()'),
                     self.do_plot_prev_pt_raw)
        self.connect(self.ui.pushButton_nextPtNumber, QtCore.SIGNAL('clicked()'),
                     self.do_plot_next_pt_raw)
        self.connect(self.ui.pushButton_showPtList, QtCore.SIGNAL('clicked()'),
                     self.show_scan_pt_list)
        self.connect(self.ui.pushButton_usePt4UB, QtCore.SIGNAL('clicked()'),
                     self.do_add_peak_to_find)
        self.connect(self.ui.pushButton_addPeakNoIndex, QtCore.SIGNAL('clicked()'),
                     self.do_add_peak_no_index)

        # Tab 'calculate ub matrix'
        self.connect(self.ui.pushButton_findPeak, QtCore.SIGNAL('clicked()'),
                     self.do_find_peak)
        self.connect(self.ui.pushButton_addPeakToCalUB, QtCore.SIGNAL('clicked()'),
                     self.do_add_ub_peak)
        self.connect(self.ui.pushButton_calUB, QtCore.SIGNAL('clicked()'),
                     self.do_cal_ub_matrix)
        self.connect(self.ui.pushButton_acceptUB, QtCore.SIGNAL('clicked()'),
                     self.doAcceptCalUB)
        self.connect(self.ui.pushButton_indexUBPeaks, QtCore.SIGNAL('clicked()'),
                     self.do_index_ub_peaks)
        self.connect(self.ui.pushButton_deleteUBPeak, QtCore.SIGNAL('clicked()'),
                     self.do_del_ub_peaks)
        self.connect(self.ui.pushButton_clearUBPeakTable, QtCore.SIGNAL('clicked()'),
                     self.do_clear_ub_peaks)
        self.connect(self.ui.pushButton_resetPeakHKLs, QtCore.SIGNAL('clicked()'),
                     self.do_reset_ub_peaks_hkl)
        self.connect(self.ui.pushButton_selectAllPeaks, QtCore.SIGNAL('clicked()'),
                     self.do_select_all_peaks)
        self.connect(self.ui.pushButton_viewScan3D, QtCore.SIGNAL('clicked()'),
                     self.do_view_data_3d)
        self.connect(self.ui.pushButton_plotSelectedData, QtCore.SIGNAL('clicked()'),
                     self.do_view_data_set_3d)

        self.connect(self.ui.pushButton_refineUB, QtCore.SIGNAL('clicked()'),
                     self.do_refine_ub)
        self.connect(self.ui.pushButton_refineUBFFT, QtCore.SIGNAL('clicked()'),
                     self.do_refine_ub_fft)

        # Tab 'Merge'
        self.connect(self.ui.pushButton_setUBSliceView, QtCore.SIGNAL('clicked()'),
                     self.do_set_ub_sv)
        self.connect(self.ui.pushButton_addScanSliceView, QtCore.SIGNAL('clicked()'),
                     self.do_add_scans_merge)
        self.connect(self.ui.pushButton_process4SliceView, QtCore.SIGNAL('clicked()'),
                     self.do_merge_scans)
        self.connect(self.ui.pushButton_readyToIntegratePeak, QtCore.SIGNAL('clicked()'),
                     self.do_advance_to_integrate_peaks)

        # Tab 'Advanced'
        self.connect(self.ui.pushButton_useDefaultDir, QtCore.SIGNAL('clicked()'),
                     self.do_setup_dir_default)
        self.connect(self.ui.pushButton_browseLocalCache, QtCore.SIGNAL('clicked()'),
                     self.do_browse_local_cache_dir)
        self.connect(self.ui.pushButton_browseWorkDir, QtCore.SIGNAL('clicked()'),
                     self.do_browse_working_dir)
        self.connect(self.ui.comboBox_instrument, QtCore.SIGNAL('currentIndexChanged(int)'),
                     self.change_instrument_name)

        # Tab 'Integrate Peaks'
        self.connect(self.ui.pushButton_integratePeak, QtCore.SIGNAL('clicked()'),
                     self.do_integrate_peaks)
        self.connect(self.ui.pushButton_findPeaks, QtCore.SIGNAL('clicked()'),
                     self.do_find_peaks_integrate)

        # Tab survey
        self.connect(self.ui.pushButton_survey, QtCore.SIGNAL('clicked()'),
                     self.do_survey)
        self.connect(self.ui.pushButton_saveSurvey, QtCore.SIGNAL('clicked()'),
                     self.do_save_survey)
        self.connect(self.ui.pushButton_loadSurvey, QtCore.SIGNAL('clicked()'),
                     self.do_load_survey)
        self.connect(self.ui.pushButton_viewSurveyPeak, QtCore.SIGNAL('clicked()'),
                     self.do_view_survey_peak)
        self.connect(self.ui.pushButton_addPeaksToRefine, QtCore.SIGNAL('clicked()'),
                     self.do_add_peaks_for_ub)

        # Menu
        self.connect(self.ui.actionExit, QtCore.SIGNAL('triggered()'),
                     self.menu_quit)

        self.connect(self.ui.actionSave_Session, QtCore.SIGNAL('triggered()'),
                     self.save_current_session)
        self.connect(self.ui.actionLoad_Session, QtCore.SIGNAL('triggered()'),
                     self.load_session)

        # Validator ... (NEXT)

        # Declaration of class variable
        # some configuration
        self._homeSrcDir = os.getcwd()
        self._homeDir = os.getcwd()

        # Control
        self._myControl = r4c.CWSCDReductionControl(self._instrument)
        self._allowDownload = True
        self._dataAccessMode = 'Download'

        # Initial setup
        self.ui.tabWidget.setCurrentIndex(0)
        self._init_table_widgets()
        self.ui.radioButton_ubFromTab1.setChecked(True)
        self.ui.lineEdit_numSurveyOutput.setText('50')
        self.ui.checkBox_loadHKLfromFile.setChecked(True)

        # Tab 'Access'
        self.ui.lineEdit_url.setText('http://neutron.ornl.gov/user_data/hb3a/')
        self.ui.comboBox_mode.setCurrentIndex(0)
        self.ui.lineEdit_localSpiceDir.setEnabled(True)
        self.ui.pushButton_browseLocalDataDir.setEnabled(True)

        # QSettings
        self.load_settings()

        return

    def closeEvent(self, QCloseEvent):
        """
        Close event
        :param QCloseEvent:
        :return:
        """
        print '[QCloseEvent=]', str(QCloseEvent)
        self.menu_quit()

    def _init_table_widgets(self):
        """ Initialize the table widgets
        :return:
        """
        # UB-peak table
        # NOTE: have to call this because pyqt set column and row to 0 after __init__
        #       thus a 2-step initialization has to been adopted
        self.ui.tableWidget_peaksCalUB.setup()
        self.ui.tableWidget_ubMatrix.setup()
        self.ui.tableWidget_ubMergeScan.setup()
        self.ui.tableWidget_surveyTable.setup()
        self.ui.tableWidget_peakIntegration.setup()
        self.ui.tableWidget_mergeScans.setup()

        return

    def do_advance_to_integrate_peaks(self):
        """
        Advance from 'merge'-tab to peak integration tab
        :return:
        """
        # Check whether there is any scan merged and selected
        ret_list = self.ui.tableWidget_mergeScans.get_rows_by_state('Done')
        if len(ret_list) == 0:
            self.pop_one_button_dialog('No scan is selected for integration!')
            return
        else:
            print '[DB] Total %d rows are selected for peak integration.' % len(ret_list)

        # Switch tab
        self.ui.tabWidget.setCurrentIndex(5)

        # Add table
        for row_index in ret_list:
            merged_ws_name = self.ui.tableWidget_mergeScans.get_merged_ws_name(row_index)
            status, merged_info = self._myControl.get_merged_scan_info(merged_ws_name)
            if status is False:
                err_msg = merged_info
                self.pop_one_button_dialog(err_msg)
                return
            else:
                print '[DB] Add selected row %d. ' % row_index
                status, msg = self.ui.tableWidget_peakIntegration.append_scan(merged_info)
                if status is False:
                    self.pop_one_button_dialog(msg)
        # END-FOR

        return

    def do_find_peaks_integrate(self):
        """ Find the centre of the 3D peak of the merged MD workspace
        :return:
        """
        row_index_list = self.ui.tableWidget_peakIntegration.get_selected_rows()

        for i_row in row_index_list:
            md_ws_name = self.ui.tableWidget_peakIntegration.get_md_ws_name(i_row)
            status, peak_centre_tuple = self._myControl.find_peak_centre_md(md_ws_name)
            if status is True:
                q_centre, hkl_centre = peak_centre_tuple
                self.ui.tableWidget_peakIntegration.set_q(i_row, q_centre)
                self.ui.tableWidget_peakIntegration.set_hkl(i_row, hkl_centre)
            # END-IF
        # END-FOR

        return

    def try_plot_3d(self):
        """ Mock test
        :return:
        """
        import numpy as np
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt

        # Dataset generation
        a, b, c = 10., 28., 8. / 3.
        def lorenz_map(X, dt = 1e-2):
            X_dt = np.array([a * (X[1] - X[0]), X[0] * (b - X[2]) - X[1], X[0] * X[1] - c * X[2]])
            return X + dt * X_dt

        def parse3DFile(filename):
            """
            """
            ifile = open(filename, 'r')
            lines = ifile.readlines()
            ifile.close()

            points = np.zeros((len(lines), 3))
            for i in xrange(len(lines)):
                line = lines[i].strip()
                terms = line.split(',')
                for j in xrange(3):
                    points[i][j] = float(terms[j])
            return points

        points = np.zeros((2000, 3))
        X = np.array([.1, .0, .0])
        for i in range(points.shape[0]):
            points[i], X = X, lorenz_map(X)

        #points0 = parse3DFile('exp355_scan38_pt11_qsample.dat')
        #points1 = parse3DFile('hkl.dat')

        # Plotting
        fig = plt.figure()
        ax = fig.gca(projection = '3d')

        ax.set_xlim(1.0, 2.1)
        ax.set_ylim(-3.2, 2.6)
        ax.set_zlim(-5.2, 2.5)

        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')
        ax.set_title('Lorenz Attractor a=%0.2f b=%0.2f c=%0.2f' % (a, b, c))

        ax.scatter(points[:, 0], points[:, 1],  points[:, 2], zdir = 'y', c = 'r')
        #ax.scatter(points0[:, 0], points0[:, 1],  points0[:, 2], zdir = 'y', c = 'r')  # c = 'k', 'r'
        #ax.scatter(points1[:, 0], points1[:, 1],  points1[:, 2], zdir = 'y', c = (0.6, 0.1, 0.9))

        plt.show()

    def do_integrate_peaks(self):
        """ Integrate peaks
        Integrate peaks
        :return:
        """
        # Determine the workspace/scan to integrate
        merged_ws_name = str(self.ui.comboBox_mergedScanWSName.currentText()).strip()
        if len(merged_ws_name) == 0:
            # merge a scan
            scan_number = gutil.parse_integers_editors([self.ui.lineEdit_scanIntegratePeak])[0]
            if scan_number is None:
                self.pop_one_button_dialog('Scan number is not given!')
                return
            pt_list = gutil.parse_integer_list(str(self.ui.lineEdit_ptNumListIntPeak.text()))
            self._myControl.merge_pts_in_scan(exp_number, scan_number, pt_list, target_ws_name=None,
                                              target_frame=None)
        else:
            if self._myControl.does_workspace_exist(merged_ws_name) is False:
                self.pop_one_button_dialog('Merged MDEventWorkspace %s does not exist!' % merged_ws_name)
                return
        # END-IF-ELSE

        """
        # Get peak integration parameters
        line_editors = [self.ui.lineEdit_peakRadius,
                        self.ui.lineEdit_bkgdInnerR,
                        self.ui.lineEdit_bkgdOuterR
                        ]
        status, value_list = gutil.parse_float_editors(line_editors)
        if status is False:
            err_msg = value_list
            print '[DB] Error message: %s' % err_msg
            return
        else:
            peak_radius = value_list[0]
            bkgd_inner_radius = value_list[1]
            bkgd_outer_radius = value_list[2]

        # Get peak integration options
        adapt_q_bkgd = self.ui.checkBox_adaptQBkgd.isChecked()
        integrate_on_edge = self.ui.checkBox_integrateOnEdge.isChecked()
        is_cylinder = self.ui.checkBox_cylinder.isChecked()

        # Choose the peaks to be integrated
        row_index_list = self.ui.tableWidget_peakIntegration.get_selected_rows()

        for i_row in row_index_list:
            md_ws_name = self.ui.tableWidget_peakIntegration.get_md_ws_name(i_row)
            exp_num = None
            scan_num = self.ui.tableWidget_peakIntegration.get_scan_number(i_row)
            pt_list = None
            self._myControl.integrate_peaks(exp_num, scan_num, pt_list, md_ws_name,
                                            peak_radius, bkgd_inner_radius, bkgd_outer_radius,
                                            is_cylinder)
        # END-FOR
        """

        return

    def do_refine_ub(self):
        """
        Refine UB matrix
        :return:
        """
        # Collecting all peaks that will be used to refine UB matrix
        row_index_list = self.ui.tableWidget_peaksCalUB.get_selected_rows(True)
        if len(row_index_list) < 3:
            err_msg = 'At least 3 peaks must be selected to refine UB matrix.' \
                      'Now it is only %d selected.' % len(row_index_list)
            self.pop_one_button_dialog(err_msg)
            return

        # loop over all peaks for peak information
        peak_info_list = list()
        status, exp_number = gutil.parse_integers_editors(self.ui.lineEdit_exp)
        assert status
        for i_row in row_index_list:
            scan_num, pt_num = self.ui.tableWidget_peaksCalUB.get_exp_info(i_row)
            try:
                peak_info = self._myControl.get_peak_info(exp_number, scan_num, pt_num)
            except AssertionError as ass_err:
                self.pop_one_button_dialog(str(ass_err))
                return
            assert isinstance(peak_info, r4c.PeakInfo)
            peak_info_list.append(peak_info)
        # END-FOR

        # Refine UB matrix
        try:
            self._myControl.refine_ub_matrix(peak_info_list)
        except AssertionError as error:
            self.pop_one_button_dialog(str(error))
            return

        # Deal with result
        ub_matrix, lattice, lattice_error = self._myControl.get_refined_ub_matrix()
        # ub matrix
        self.ui.tableWidget_ubMatrix.set_from_matrix(ub_matrix)

        # lattice parameter
        assert isinstance(lattice, list)
        assert len(lattice) == 6
        self.ui.lineEdit_aUnitCell.setText('%.5f' % lattice[0])
        self.ui.lineEdit_bUnitCell.setText('%.5f' % lattice[1])
        self.ui.lineEdit_cUnitCell.setText('%.5f' % lattice[2])
        self.ui.lineEdit_alphaUnitCell.setText('%.5f' % lattice[3])
        self.ui.lineEdit_betaUnitCell.setText('%.5f' % lattice[4])
        self.ui.lineEdit_gammaUnitCell.setText('%.5f' % lattice[5])

        assert isinstance(lattice_error, list)
        assert len(lattice_error) == 6
        self.ui.lineEdit_aError.setText('%.5f' % lattice_error[0])
        self.ui.lineEdit_bError.setText('%.5f' % lattice_error[1])
        self.ui.lineEdit_cError.setText('%.5f' % lattice_error[2])
        self.ui.lineEdit_alphaError.setText('%.5f' % lattice_error[3])
        self.ui.lineEdit_betaError.setText('%.5f' % lattice_error[4])
        self.ui.lineEdit_gammaError.setText('%.5f' % lattice_error[5])

        # TODO/NOW/1st: need to offer users with different types of UB matrix refinement tool!
        # call mantid.FindUBUsingIndexedPeaks()
        # refer to Calculate UB matrix to build PeakWorkspace

        raise RuntimeError('Next Release')

    def do_refine_ub_fft(self):
        """
        Refine UB matrix by calling FFT method
        :return:
        """
        # TODO/NOW/Combine with do_refine_ub()

        return

    def change_data_access_mode(self):
        """ Change data access mode between downloading from server and local
        Event handling methods
        :return:
        """
        new_mode = str(self.ui.comboBox_mode.currentText())
        self._dataAccessMode = new_mode

        if new_mode.startswith('Local') is True:
            self.ui.lineEdit_localSpiceDir.setEnabled(True)
            self.ui.pushButton_browseLocalDataDir.setEnabled(True)
            self.ui.lineEdit_url.setEnabled(False)
            self.ui.lineEdit_localSrcDir.setEnabled(False)
            self.ui.pushButton_browseLocalCache.setEnabled(False)
            self._allowDownload = False
        else:
            self.ui.lineEdit_localSpiceDir.setEnabled(False)
            self.ui.pushButton_browseLocalDataDir.setEnabled(False)
            self.ui.lineEdit_url.setEnabled(True)
            self.ui.lineEdit_localSrcDir.setEnabled(True)
            self.ui.pushButton_browseLocalCache.setEnabled(True)
            self._allowDownload = True

        return

    def change_instrument_name(self):
        """ Handing the event as the instrument name is changed
        :return:
        """
        new_instrument = str(self.ui.comboBox_instrument.currentText())
        self.pop_one_button_dialog('Change of instrument during data processing is dangerous.')
        status, error_message = self._myControl.set_instrument_name(new_instrument)
        if status is False:
            self.pop_one_button_dialog(error_message)

        return

    def do_add_peak_no_index(self):
        """
        Purpose: add a peak from 'View Raw Data' tab to the UB peak table
        without indexing it
        :return:
        """
        # Get exp, scan and Pt information
        status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_exp,
                                                        self.ui.lineEdit_run,
                                                        self.ui.lineEdit_rawDataPtNo])
        if not status:
            err_msg = ret_obj
            self.pop_one_button_dialog(err_msg)
            return
        else:
            int_list = ret_obj
            exp_no, scan_no, pt_no = int_list

        # Switch tab
        self.ui.tabWidget.setCurrentIndex(2)

        # Find and index peak
        self._myControl.find_peak(exp_no, scan_no, [pt_no])
        try:
            peak_info = self._myControl.get_peak_info(exp_no, scan_no, pt_no)
            self.set_ub_peak_table(peak_info)
        except AssertionError as ass_err:
            self.pop_one_button_dialog(str(ass_err))
            return

        return

    def do_add_ub_peak(self):
        """ Add current to ub peaks
        :return:
        """
        # Add peak
        status, int_list = gutil.parse_integers_editors([self.ui.lineEdit_exp,
                                                         self.ui.lineEdit_scanNumber])
        if status is False:
            self.pop_one_button_dialog(int_list)
            return
        exp_no, scan_no = int_list

        # Get HKL from GUI
        status, float_list = gutil.parse_float_editors([self.ui.lineEdit_H,
                                                        self.ui.lineEdit_K,
                                                        self.ui.lineEdit_L])
        if status is False:
            err_msg = float_list
            self.pop_one_button_dialog(err_msg)
            return
        h, k, l = float_list

        try:
            peak_info_obj = self._myControl.get_peak_info(exp_no, scan_no)
        except AssertionError as ass_err:
            self.pop_one_button_dialog(str(ass_err))
            return

        assert isinstance(peak_info_obj, r4c.PeakInfo)
        if self.ui.checkBox_roundHKLInt.isChecked():
            h = math.copysign(1, h)*int(abs(h)+0.5)
            k = math.copysign(1, k)*int(abs(k)+0.5)
            l = math.copysign(1, l)*int(abs(l)+0.5)
        peak_info_obj.set_user_hkl(h, k, l)
        self.set_ub_peak_table(peak_info_obj)

        # Clear
        self.ui.lineEdit_scanNumber.setText('')

        self.ui.lineEdit_sampleQx.setText('')
        self.ui.lineEdit_sampleQy.setText('')
        self.ui.lineEdit_sampleQz.setText('')

        self.ui.lineEdit_H.setText('')
        self.ui.lineEdit_K.setText('')
        self.ui.lineEdit_L.setText('')

        return

    def doAcceptCalUB(self):
        """ Accept the calculated UB matrix
        """
        raise RuntimeError('ASAP')

    def doAddScanPtToRefineUB(self):
        """ Add scan/pt numbers to the list of data points for refining ub matrix

        And the added scan number and pt numbers will be reflected in the (left sidebar)

        """
        raise RuntimeError("ASAP")

    def do_add_peak_to_find(self):
        """
        Add the scan/pt to the next
        :return:
        """
        # peak finding will use all points in the selected scan.
        scan_no = self.ui.lineEdit_run.text()
        self.ui.tabWidget.setCurrentIndex(2)
        self.ui.lineEdit_scanNumber.setText(scan_no)

        return

    def do_add_peaks_for_ub(self):
        """ In tab-survey, merge selected scans, find peaks in merged data and
         switch to UB matrix calculation tab and add to table
        :return:
        """
        # get selected scans
        selected_row_index_list = self.ui.tableWidget_surveyTable.get_selected_rows(True)
        scan_number_list = self.ui.tableWidget_surveyTable.get_scan_numbers(selected_row_index_list)
        if len(scan_number_list) == 0:
            self.pop_one_button_dialog('No scan is selected.')
            return

        # get experiment number
        status, exp_number = gutil.parse_integers_editors(self.ui.lineEdit_exp)
        assert status

        # switch to tab-3
        self.ui.tabWidget.setCurrentIndex(2)

        # find peak and add peak
        for scan_number in scan_number_list:
            # merge peak and find peak
            self._myControl.merge_pts_in_scan(exp_number, scan_number, [], 'q-sample')
            self._myControl.find_peak(exp_number, scan_number)

            # get PeakInfo
            peak_info = self._myControl.get_peak_info(exp_number, scan_number)
            assert isinstance(peak_info, r4c.PeakInfo)

            # retrieve and set HKL from spice table
            peak_info.retrieve_hkl_from_spice_table()

            # add to table
            self.set_ub_peak_table(peak_info)
        # END-FOR

        return

    def do_browse_local_cache_dir(self):
        """ Browse local cache directory
        :return:
        """
        local_cache_dir = str(QtGui.QFileDialog.getExistingDirectory(self,
                                                                     'Get Local Cache Directory',
                                                                     self._homeSrcDir))

        # Set local directory to control
        status, error_message = self._myControl.set_local_data_dir(local_cache_dir)
        if status is False:
            self.pop_one_button_dialog(error_message)
            return

        # Synchronize to local data/spice directory and local cache directory
        if str(self.ui.lineEdit_localSpiceDir.text()) != '':
            prev_dir = str(self.ui.lineEdit_localSrcDir.text())
            self.pop_one_button_dialog('Local data directory was set up as %s' %
                                       prev_dir)
        self.ui.lineEdit_localSrcDir.setText(local_cache_dir)
        self.ui.lineEdit_localSpiceDir.setText(local_cache_dir)

        return

    def do_browse_local_spice_data(self):
        """ Browse local source SPICE data directory
        """
        src_spice_dir = str(QtGui.QFileDialog.getExistingDirectory(self, 'Get Directory',
                                                                   self._homeSrcDir))
        # Set local data directory to controller
        status, error_message = self._myControl.set_local_data_dir(src_spice_dir)
        if status is False:
            self.pop_one_button_dialog(error_message)
            return

        self._homeSrcDir = src_spice_dir
        self.ui.lineEdit_localSpiceDir.setText(src_spice_dir)

        return

    def do_browse_working_dir(self):
        """
        Browse and set up working directory
        :return:
        """
        work_dir = str(QtGui.QFileDialog.getExistingDirectory(self, 'Get Working Directory', self._homeDir))
        status, error_message = self._myControl.set_working_directory(work_dir)
        if status is False:
            self.pop_one_button_dialog(error_message)
        else:
            self.ui.lineEdit_workDir.setText(work_dir)

        return

    def do_cal_ub_matrix(self):
        """ Calculate UB matrix by 2 or 3 reflections
        """
        # Get reflections selected to calculate UB matrix
        num_rows = self.ui.tableWidget_peaksCalUB.rowCount()
        peak_info_list = list()
        status, exp_number = gutil.parse_integers_editors(self.ui.lineEdit_exp)
        assert status
        for i_row in xrange(num_rows):
            if self.ui.tableWidget_peaksCalUB.is_selected(i_row) is True:
                scan_num, pt_num = self.ui.tableWidget_peaksCalUB.get_exp_info(i_row)
                if pt_num < 0:
                    pt_num = None
                peak_info = self._myControl.get_peak_info(exp_number, scan_num, pt_num)
                assert isinstance(peak_info, r4c.PeakInfo)
                peak_info_list.append(peak_info)
        # END-FOR

        # Get lattice
        status, ret_obj = self._get_lattice_parameters()
        if status is True:
            a, b, c, alpha, beta, gamma = ret_obj
        else:
            err_msg = ret_obj
            self.pop_one_button_dialog(err_msg)
            return

        # Calculate UB matrix
        status, ub_matrix = self._myControl.calculate_ub_matrix(peak_info_list, a, b, c,
                                                                alpha, beta, gamma)

        # Deal with result
        if status is True:
            self.ui.tableWidget_ubMatrix.set_from_matrix(ub_matrix)

        else:
            err_msg = ub_matrix
            self.pop_one_button_dialog(err_msg)

        return

    def do_clear_ub_peaks(self):
        """
        Clear all peaks in UB-Peak table
        :return:
        """
        num_rows = self.ui.tableWidget_peaksCalUB.rowCount()
        row_number_list = range(num_rows)
        self.ui.tableWidget_peaksCalUB.delete_rows(row_number_list)

        return

    def do_del_ub_peaks(self):
        """
        Delete a peak in UB-Peak table
        :return:
        """
        # Find out the lines to get deleted
        row_num_list = self.ui.tableWidget_peaksCalUB.get_selected_rows()
        print '[DB] Row %s are selected' % str(row_num_list)

        # Delete
        self.ui.tableWidget_peaksCalUB.delete_rows(row_num_list)

        return

    def do_download_spice_data(self):
        """ Download SPICE data
        :return:
        """
        # Check scans to download
        scan_list_str = str(self.ui.lineEdit_downloadScans.text())
        if len(scan_list_str) > 0:
            # user specifies scans to download
            valid, scan_list = fcutil.parse_int_array(scan_list_str)
            if valid is False:
                error_message = scan_list
                self.pop_one_button_dialog(error_message)
        else:
            # Get all scans
            status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_exp])
            if status is False:
                self.pop_one_button_dialog(ret_obj)
                return
            exp_no = ret_obj
            assert isinstance(exp_no, int)
            server_url = str(self.ui.lineEdit_url.text())
            scan_list = fcutil.get_scans_list(server_url, exp_no, return_list=True)
        self.pop_one_button_dialog('Going to download scans %s.' % str(scan_list))

        # Check location
        destination_dir = str(self.ui.lineEdit_localSrcDir.text())
        status, error_message = self._myControl.set_local_data_dir(destination_dir)
        if status is False:
            self.pop_one_button_dialog(error_message)
        else:
            self.pop_one_button_dialog('Spice files will be downloaded to %s.' % destination_dir)

        # Set up myControl for downloading data
        exp_no = int(self.ui.lineEdit_exp.text())
        self._myControl.set_exp_number(exp_no)

        server_url = str(self.ui.lineEdit_url.text())
        status, error_message = self._myControl.set_server_url(server_url)
        if status is False:
            self.pop_one_button_dialog(error_message)
            return

        # Download
        self._myControl.download_data_set(scan_list)

        return

    def do_find_peak(self):
        """ Find peak in a given scan and record it
        """
        # Get experiment, scan and pt
        status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_exp,
                                                        self.ui.lineEdit_scanNumber])
        if status is True:
            exp_no, scan_no = ret_obj
        else:
            self.pop_one_button_dialog(ret_obj)
            return

        # merge peak if necessary
        if self._myControl.has_merged_data(exp_no, scan_no) is False:
            status, err_msg = self._myControl.merge_pts_in_scan(exp_no, scan_no, [], 'q-sample')
            if status is False:
                self.pop_one_button_dialog(err_msg)

        # Find peak
        status, err_msg = self._myControl.find_peak(exp_no, scan_no)
        if status is False:
            self.pop_one_button_dialog(ret_obj)
            return

        # Get information from the latest (integrated) peak
        if self.ui.checkBox_loadHKLfromFile.isChecked() is True:
            # This is the first time that in the workflow to get HKL from MD workspace
            peak_info = self._myControl.get_peak_info(exp_no, scan_no)
            try:
                peak_info.retrieve_hkl_from_spice_table()
            except RuntimeError as run_err:
                self.pop_one_button_dialog('Unable to locate peak info due to %s.' % str(run_err))
        # END-IF

        # Set up correct values to table tableWidget_peaksCalUB
        peak_info = self._myControl.get_peak_info(exp_no, scan_no)
        h, k, l = peak_info.get_user_hkl()
        self.ui.lineEdit_H.setText('%.2f' % h)
        self.ui.lineEdit_K.setText('%.2f' % k)
        self.ui.lineEdit_L.setText('%.2f' % l)

        q_x, q_y, q_z = peak_info.get_peak_centre()
        self.ui.lineEdit_sampleQx.setText('%.5E' % q_x)
        self.ui.lineEdit_sampleQy.setText('%.5E' % q_y)
        self.ui.lineEdit_sampleQz.setText('%.5E' % q_z)

        return

    def do_index_ub_peaks(self):
        """ Index the peaks in the UB matrix peak table
        :return:
        """
        # Get UB matrix
        ub_matrix = self.ui.tableWidget_ubMatrix.get_matrix()
        print '[Info] Get UB matrix from table ', ub_matrix

        # Index all peaks
        num_peaks = self.ui.tableWidget_peaksCalUB.rowCount()
        err_msg = ''
        for i_peak in xrange(num_peaks):
            scan_no, pt_no = self.ui.tableWidget_peaksCalUB.get_exp_info(i_peak)
            status, ret_obj = self._myControl.index_peak(ub_matrix, scan_number=scan_no)
            if status is True:
                hkl_value = ret_obj[0]
                hkl_error = ret_obj[1]
                self.ui.tableWidget_peaksCalUB.set_hkl(i_peak, hkl_value, hkl_error)
            else:
                err_msg += ret_obj + '\n'
        # END-FOR

        if len(err_msg) > 0:
            self.pop_one_button_dialog(err_msg)

        return

    def do_list_scans(self):
        """ List all scans available
        :return:
        """
        # Experiment number
        exp_no = int(self.ui.lineEdit_exp.text())

        access_mode = str(self.ui.comboBox_mode.currentText())
        if access_mode == 'Local':
            spice_dir = str(self.ui.lineEdit_localSpiceDir.text())
            message = fcutil.get_scans_list_local_disk(spice_dir, exp_no)
        else:
            url = str(self.ui.lineEdit_url.text())
            message = fcutil.get_scans_list(url, exp_no)

        self.pop_one_button_dialog(message)

        return

    def do_load_scan_info(self):
        """ Load SIICE's scan file
        :return:
        """
        # Get scan number
        status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_run])
        if status is True:
            scan_no = ret_obj[0]
        else:
            err_msg = ret_obj
            self.pop_one_button_dialog('Unable to get scan number in raw data tab due to %s.' % err_msg)
            return

        status, err_msg = self._myControl.load_spice_scan_file(exp_no=None, scan_no=scan_no)
        if status is False:
            self.pop_one_button_dialog(err_msg)

        return

    def do_load_survey(self):
        """ Load csv file containing experiment-scan survey's result.
        :return:
        """
        # check validity
        num_rows = int(self.ui.lineEdit_numSurveyOutput.text())

        # get the csv file
        file_filter = 'CSV Files (*.csv);;All Files (*.*)'
        csv_file_name = str(QtGui.QFileDialog.getOpenFileName(self, 'Open Exp-Scan Survey File', self._homeDir,
                                                              file_filter))

        # call controller to load
        survey_tuple = self._myControl.load_scan_survey_file(csv_file_name)
        scan_sum_list = survey_tuple[1]
        assert isinstance(scan_sum_list, list), 'Returned value from load scan survey file must be a dictionary.'

        # set the table
        self.ui.tableWidget_surveyTable.set_survey_result(scan_sum_list)
        self.ui.tableWidget_surveyTable.remove_all_rows()
        self.ui.tableWidget_surveyTable.show_reflections(num_rows)

        return

    def do_plot_pt_raw(self):
        """ Plot the Pt.
        """
        # Get measurement pt and the file number
        status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_exp,
                                                        self.ui.lineEdit_run,
                                                        self.ui.lineEdit_rawDataPtNo])
        if status is True:
            exp_no = ret_obj[0]
            scan_no = ret_obj[1]
            pt_no = ret_obj[2]
        else:
            self.pop_one_button_dialog(ret_obj)
            return

        # Call to plot 2D
        self._plot_raw_xml_2d(exp_no, scan_no, pt_no)

        # Information
        info = '%-10s: %d\n%-10s: %d\n%-10s: %d\n' % ('Exp', exp_no,
                                                      'Scan', scan_no,
                                                      'Pt', pt_no)
        self.ui.plainTextEdit_rawDataInformation.setPlainText(info)

        return

    def do_plot_prev_pt_raw(self):
        """ Plot the Pt.
        """
        # Get measurement pt and the file number
        status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_exp,
                                                        self.ui.lineEdit_run,
                                                        self.ui.lineEdit_rawDataPtNo])
        if status is True:
            exp_no = ret_obj[0]
            scan_no = ret_obj[1]
            pt_no = ret_obj[2]
        else:
            self.pop_one_button_dialog(ret_obj)
            return

        # Previous one
        pt_no -= 1
        if pt_no <= 0:
            self.pop_one_button_dialog('Pt. = 1 is the first one.')
            return
        else:
            self.ui.lineEdit_rawDataPtNo.setText('%d' % pt_no)

        # Plot
        self._plot_raw_xml_2d(exp_no, scan_no, pt_no)

        return

    def do_plot_next_pt_raw(self):
        """ Plot the Pt.
        """
        # Get measurement pt and the file number
        status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_exp,
                                                        self.ui.lineEdit_run,
                                                        self.ui.lineEdit_rawDataPtNo])
        if status is True:
            exp_no = ret_obj[0]
            scan_no = ret_obj[1]
            pt_no = ret_obj[2]
        else:
            self.pop_one_button_dialog(ret_obj)
            return

        # Previous one
        pt_no += 1
        # get last Pt. number
        status, last_pt_no = self._myControl.get_pt_numbers(exp_no, scan_no)
        if status is False:
            error_message = last_pt_no
            self.pop_one_button_dialog('Unable to access Spice table for scan %d. Reason" %s.' % (
                scan_no, error_message))
        if pt_no > last_pt_no:
            self.pop_one_button_dialog('Pt. = %d is the last one of scan %d.' % (pt_no, scan_no))
            return
        else:
            self.ui.lineEdit_rawDataPtNo.setText('%d' % pt_no)

        # Plot
        self._plot_raw_xml_2d(exp_no, scan_no, pt_no)

        return

    def do_add_scans_merge(self):
        """ Add scans to merge
        :return:
        """
        # Get list of scans
        scan_list = gutil.parse_integer_list(str(self.ui.lineEdit_listScansSliceView.text()))
        if len(scan_list) == 0:
            self.pop_one_button_dialog('Scan list is empty.')

        # Set table
        self.ui.tableWidget_mergeScans.append_scans(scans=scan_list)

        return

    def do_merge_scans(self):
        """ Process data for slicing view
        :return:
        """
        # Get UB matrix
        ub_matrix = self.ui.tableWidget_ubMergeScan.get_matrix()
        self._myControl.set_ub_matrix(exp_number=None, ub_matrix=ub_matrix)

        # Warning
        self.pop_one_button_dialog('Data processing is long. Be patient!')

        # Process
        base_name = str(self.ui.lineEdit_baseMergeMDName.text())
        scan_row_list = self.ui.tableWidget_mergeScans.get_scan_list()
        print '[DB] %d scans have been selected to merge.' % len(scan_row_list)
        frame = str(self.ui.comboBox_mergeScanFrame.currentText())
        for tup2 in scan_row_list:
            #
            scan_no, i_row = tup2

            # Download/check SPICE file
            self._myControl.download_spice_file(None, scan_no, over_write=False)

            # Get some information
            status, pt_list = self._myControl.get_pt_numbers(None, scan_no, load_spice_scan=True)
            if status is False:
                err_msg = pt_list
                self.pop_one_button_dialog('Failed to get Pt. number: %s' % err_msg)
                return
            else:
                # Set information to table
                err_msg = self.ui.tableWidget_mergeScans.set_scan_pt(scan_no, pt_list)
                if len(err_msg) > 0:
                    self.pop_one_button_dialog(err_msg)

            self.ui.tableWidget_mergeScans.set_status_by_row(i_row, 'In Processing')
            merge_status = 'UNKNOWN'
            merged_name = '???'
            group_name = '???'

            status, ret_tup = self._myControl.merge_pts_in_scan(exp_no=None, scan_no=scan_no,
                                                                pt_num_list=[], target_frame=frame)
            merge_status = 'Done'
            merged_name = ret_tup[0]
            group_name = ret_tup[1]

            if status is False:
                merge_status = 'Failed. Reason: %s' % str(e)
                merged_name = ''
                group_name = ''
                print merge_status
            else:
                self.ui.tableWidget_mergeScans.set_status_by_row(i_row, merge_status)
                self.ui.tableWidget_mergeScans.set_ws_names_by_row(i_row, merged_name, group_name)

            # Sleep for a while
            time.sleep(0.1)
        # END-FOR

        return

    def do_reset_ub_peaks_hkl(self):
        """
        Reset user specified HKL value to peak table
        :return:
        """
        num_rows = self.ui.tableWidget_peaksCalUB.rowCount()
        for i_row in xrange(num_rows):
            print '[DB] Update row %d' % (i_row)
            scan, pt = self.ui.tableWidget_peaksCalUB.get_scan_pt(i_row)
            peak_info = self._myControl.get_peak_info(None, scan, pt)
            h, k, l = peak_info.get_user_hkl()
            self.ui.tableWidget_peaksCalUB.update_hkl(i_row, h, k, l)
        # END-FOR

        return

    def do_save_survey(self):
        """
        Save the survey to a file
        :return:
        """
        # Get file name
        file_filter = 'CSV Files (*.csv);;All Files (*.*)'
        out_file_name = str(QtGui.QFileDialog.getSaveFileName(self, 'Save scan survey result',
                                                              self._homeDir, file_filter))

        # Save file
        self._myControl.save_scan_survey(out_file_name)

        return

    def do_select_all_peaks(self):
        """
        Purpose: select all peaks in table tableWidget_peaksCalUB
        :return:
        """
        # TODO/NOW/1st: Implement ASAP

        return

    def do_set_experiment(self):
        """ Set experiment
        :return:
        """
        status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_exp])
        if status is True:
            exp_number = ret_obj[0]
            curr_exp_number = self._myControl.get_experiment()
            if curr_exp_number is not None and exp_number != curr_exp_number:
                self.pop_one_button_dialog('Changing experiment to %d.  Clean previous experiment %d\'s result'
                                           ' in Mantid manually.' % (exp_number, curr_exp_number))
            self._myControl.set_exp_number(exp_number)
            self.ui.lineEdit_exp.setStyleSheet('color: black')
        else:
            err_msg = ret_obj
            self.pop_one_button_dialog('Unable to set experiment as %s' % err_msg)
            self.ui.lineEdit_exp.setStyleSheet('color: red')

        self.ui.tabWidget.setCurrentIndex(0)

        return

    def do_set_ub_sv(self):
        """ Set UB matrix in Slice view
        :return:
        """
        if self.ui.radioButton_ubFromTab1.isChecked():
            # from tab 'calculate UB'
            self.ui.tableWidget_ubMergeScan.set_from_matrix(self.ui.tableWidget_ubMatrix.get_matrix())
        elif self.ui.radioButton_ubFromTab3.isChecked():
            # from tab 'refine UB'
            self.ui.tableWidget_ubMergeScan.set_from_matrix(self.ui.tableWidget_refinedUB.get_matrix())
        elif self.ui.radioButton_ubFromList.isChecked():
            # set ub matrix manually
            ub_str = str(self.ui.plainTextEdit_ubInput.toPlainText())
            status, ret_obj = gutil.parse_float_array(ub_str)
            if status is False:
                # unable to parse to float arrays
                self.pop_one_button_dialog(ret_obj)
            elif len(ret_obj) != 9:
                # number of floats is not 9
                self.pop_one_button_dialog('Requiring 9 floats for UB matrix.  Only %d are given.' % len(ret_obj))
            else:
                # in good UB matrix format
                ub_str = ret_obj
                option = str(self.ui.comboBox_ubOption.currentText())
                if option.lower().count('spice') > 0:
                    # convert from SPICE UB matrix to Mantid one
                    spice_ub = gutil.convert_str_to_matrix(ub_str, (3, 3))
                    mantid_ub = r4c.convert_spice_ub_to_mantid(spice_ub)
                    self.ui.tableWidget_ubMergeScan.set_from_matrix(mantid_ub)
                else:
                    self.ui.tableWidget_ubMergeScan.set_from_list(ub_str)
        else:
            self.pop_one_button_dialog('None is selected to set UB matrix.')

        return

    def do_setup_dir_default(self):
        """
        Set up default directory for storing data and working
        :return:
        """
        home_dir = os.path.expanduser('~')

        # Data cache directory
        data_cache_dir = os.path.join(home_dir, 'Temp/HB3ATest')
        self.ui.lineEdit_localSpiceDir.setText(data_cache_dir)
        self.ui.lineEdit_localSrcDir.setText(data_cache_dir)

        work_dir = os.path.join(data_cache_dir, 'Workspace')
        self.ui.lineEdit_workDir.setText(work_dir)

        return

    def do_apply_setup(self):
        """
        Apply set up ...
        :return:
        """
        # Local data directory
        local_data_dir = str(self.ui.lineEdit_localSpiceDir.text())
        if os.path.exists(local_data_dir) is False:
            try:
                os.mkdir(local_data_dir)
            except OSError as os_error:
                self.pop_one_button_dialog('Unable to create local data directory %s due to %s.' % (
                    local_data_dir, str(os_error)))
                self.ui.lineEdit_localSpiceDir.setStyleSheet("color: red;")
                return
            else:
                self.ui.lineEdit_localSpiceDir.setStyleSheet("color: black;")
        # END-IF

        # Working directory
        working_dir = str(self.ui.lineEdit_workDir.text())
        if os.path.exists(working_dir) is False:
            try:
                os.mkdir(working_dir)
            except OSError as os_error:
                self.pop_one_button_dialog('Unable to create working directory %s due to %s.' % (
                    working_dir, str(os_error)))
                self.ui.lineEdit_workDir.setStyleSheet("color: red;")
                return
            else:
                self.ui.lineEdit_workDir.setStyleSheet("color: black;")
        # END-IF

        # Server URL
        data_server = str(self.ui.lineEdit_url.text())
        url_is_good = self.do_test_url()
        if url_is_good is False:
            self.ui.lineEdit_url.setStyleSheet("color: red;")
            return
        else:
            self.ui.lineEdit_url.setStyleSheet("color: black;")

        # Set to control
        self._myControl.set_local_data_dir(local_data_dir)
        self._myControl.set_working_directory(working_dir)
        self._myControl.set_server_url(data_server)

        return

    def do_survey(self):
        """
        Purpose: survey for the strongest reflections
        :return:
        """
        # Get experiment number
        exp_number = int(self.ui.lineEdit_exp.text())
        status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_surveyStartPt,
                                                        self.ui.lineEdit_surveyEndPt])
        if status is False:
            err_msg = ret_obj
            self.pop_one_button_dialog(err_msg)
        start_scan = ret_obj[0]
        end_scan = ret_obj[1]

        max_number = int(self.ui.lineEdit_numSurveyOutput.text())

        # Get value
        status, ret_obj = self._myControl.survey(exp_number, start_scan, end_scan)
        if status is False:
            self.pop_one_button_dialog(ret_obj)
            return
        scan_sum_list = ret_obj
        self.ui.tableWidget_surveyTable.set_survey_result(scan_sum_list)
        self.ui.tableWidget_surveyTable.show_reflections(max_number)

        return

    def do_test_url(self):
        """ Test whether the root URL provided specified is good
        """
        url = str(self.ui.lineEdit_url.text())

        url_is_good, err_msg = fcutil.check_url(url)
        if url_is_good is True:
            self.pop_one_button_dialog("URL %s is valid." % url)
        else:
            self.pop_one_button_dialog(err_msg)

        return url_is_good

    def do_view_data_set_3d(self):
        """

        :return:
        """
        # TODO/NOW/ Think of share codes with do_view_data_3d()
        import plot3dwindow

        self.plot3dform = plot3dwindow.Plot3DWindow(self)
        self.plot3dform.show()

        return

    def do_view_data_3d(self):
        """
        View merged scan data in 3D after FindPeaks
        :return:
        """
        # get experiment and scan number
        status, ret_obj = gutil.parse_integers_editors([self.ui.lineEdit_exp, self.ui.lineEdit_scanNumber])
        if status:
            exp_number = ret_obj[0]
            scan_number = ret_obj[1]
        else:
            self.pop_one_button_dialog(ret_obj)
            return

        # Check
        if self._myControl.has_merged_data(exp_number, scan_number) is False:
            self.pop_one_button_dialog('No merged MD workspace found for %d, %d' % (exp_number, scan_number))
            return

        # Generate data by writing out temp file
        base_file_name = 'md_%d.dat' % random.randint(1, 1001)
        md_file_name = self._myControl.export_md_data(exp_number, scan_number, base_file_name)
        peak_info = self._myControl.get_peak_info(exp_number, scan_number)
        weight_peak_centers = peak_info.get_weighted_peak_centres()
        avg_peak_centre = peak_info.get_peak_centre()

        print 'Write file to %s' % md_file_name
        for i_peak in xrange(len(weight_peak_centers)):
            peak_i = weight_peak_centers[i_peak]
            print '%f, %f, %f, %f' % (peak_i[0], peak_i[1], peak_i[2], peak_i[3])
        print
        print avg_peak_centre

        # Plot

        raise

    def do_view_survey_peak(self):
        """ View selected peaks from survey table
        Requirements: one and only 1 run is selected
        Guarantees: the scan number and pt number that are selected will be set to
            tab 'View Raw' and the tab is switched.
        :return:
        """
        # get values
        try:
            scan_num, pt_num = self.ui.tableWidget_surveyTable.get_selected_run_surveyed()
        except RuntimeError, err:
            self.pop_one_button_dialog(str(err))
            return

        # clear selection
        self.ui.tableWidget_surveyTable.select_all_rows(False)

        # switch tab
        self.ui.tabWidget.setCurrentIndex(1)
        self.ui.lineEdit_run.setText(str(scan_num))
        self.ui.lineEdit_rawDataPtNo.setText(str(pt_num))

        return

    def pop_one_button_dialog(self, message):
        """ Pop up a one-button dialog
        :param message:
        :return:
        """
        assert isinstance(message, str)
        QtGui.QMessageBox.information(self, '4-circle Data Reduction', message)

        return

    def save_current_session(self, filename=None):
        """ Save current session/value setup to
        :return:
        """
        # Set up dictionary
        save_dict = dict()

        # Setup
        save_dict['lineEdit_localSpiceDir'] = str(self.ui.lineEdit_localSpiceDir.text())
        save_dict['lineEdit_url'] = str(self.ui.lineEdit_url.text())
        save_dict['lineEdit_workDir']= str(self.ui.lineEdit_workDir.text())

        # Experiment
        save_dict['lineEdit_exp'] = str(self.ui.lineEdit_exp.text())
        save_dict['lineEdit_scanNumber'] = self.ui.lineEdit_scanNumber.text()
        save_dict['lineEdit_ptNumber'] = str(self.ui.lineEdit_ptNumber.text())

        # Lattice
        save_dict['lineEdit_a'] = str(self.ui.lineEdit_a.text())
        save_dict['lineEdit_b'] = str(self.ui.lineEdit_b.text())
        save_dict['lineEdit_c'] = str(self.ui.lineEdit_c.text())
        save_dict['lineEdit_alpha'] = str(self.ui.lineEdit_alpha.text())
        save_dict['lineEdit_beta'] = str(self.ui.lineEdit_beta.text())
        save_dict['lineEdit_gamma'] = str(self.ui.lineEdit_gamma.text())

        # Merge scan
        save_dict['plainTextEdit_ubInput'] = str(self.ui.plainTextEdit_ubInput.toPlainText())
        save_dict['lineEdit_listScansSliceView'] = str(self.ui.lineEdit_listScansSliceView.text())
        save_dict['lineEdit_baseMergeMDName'] = str(self.ui.lineEdit_baseMergeMDName.text())

        # Save to csv file
        if filename is None:
            filename = 'session_backup.csv'
        ofile = open(filename, 'w')
        writer = csv.writer(ofile)
        for key, value in save_dict.items():
            writer.writerow([key, value])
        ofile.close()

        return

    def load_session(self, filename=None):
        """
        To load a session, i.e., read it back:
        :param filename:
        :return:
        """
        if filename is None:
            filename = 'session_backup.csv'

        in_file = open(filename, 'r')
        reader = csv.reader(in_file)
        my_dict = dict(x for x in reader)

        # ...
        for key, value in my_dict.items():
            if key.startswith('lineEdit') is True:
                self.ui.__getattribute__(key).setText(value)
            elif key.startswith('plainText') is True:
                self.ui.__getattribute__(key).setPlainText(value)
            elif key.startswith('comboBox') is True:
                self.ui.__getattribute__(key).setCurrentIndex(int(value))
            else:
                self.pop_one_button_dialog('Error! Widget name %s is not supported' % key)
        # END-FOR

        # ...
        self._myControl.set_local_data_dir(str(self.ui.lineEdit_localSpiceDir.text()))

        return

    def menu_quit(self):
        """

        :return:
        """
        self.save_settings()
        self.close()

    def show_scan_pt_list(self):
        """ Show the range of Pt. in a scan
        :return:
        """
        # Get parameters
        status, inp_list = gutil.parse_integers_editors([self.ui.lineEdit_exp, self.ui.lineEdit_run])
        if status is False:
            self.pop_one_button_dialog(inp_list)
            return
        else:
            exp_no = inp_list[0]
            scan_no = inp_list[1]

        status, ret_obj = self._myControl.get_pt_numbers(exp_no, scan_no)

        # Form message
        if status is False:
            # Failed to get Pt. list
            error_message = ret_obj
            self.pop_one_button_dialog(error_message)
        else:
            # Form message
            pt_list = sorted(ret_obj)
            num_pts = len(pt_list)
            info = 'Exp %d Scan %d has %d Pt. ranging from %d to %d.\n' % (exp_no, scan_no, num_pts,
                                                                           pt_list[0], pt_list[-1])
            num_miss_pt = pt_list[-1] - pt_list[0] + 1 - num_pts
            if num_miss_pt > 0:
                info += 'There are %d Pt. skipped.\n' % num_miss_pt

            self.pop_one_button_dialog(info)

        return

    def set_ub_peak_table(self, peakinfo):
        """
        TODO/NOW/DOC
        :param peak_info:
        :return:
        """
        assert isinstance(peakinfo, r4c.PeakInfo)

        # Get data
        exp_number, scan_number = peakinfo.getExpInfo()
        h, k, l = peakinfo.get_user_hkl()
        q_x, q_y, q_z = peakinfo.get_peak_centre()
        m1 = self._myControl.get_sample_log_value(exp_number, scan_number, 1, '_m1')

        # Set to table
        status, err_msg = self.ui.tableWidget_peaksCalUB.append_row(
            [scan_number, -1, h, k, l, q_x, q_y, q_z, False, m1, ''])
        if status is False:
            self.pop_one_button_dialog(err_msg)

        return

    def save_settings(self):
        """
        Save settings (parameter set) upon quiting
        :return:
        """
        settings = QtCore.QSettings()

        # directories
        local_spice_dir = str(self.ui.lineEdit_localSpiceDir.text())
        settings.setValue("local_spice_dir", local_spice_dir)
        work_dir = str(self.ui.lineEdit_workDir.text())
        settings.setValue('work_dir', work_dir)

        # experiment number
        exp_num = str(self.ui.lineEdit_exp.text())
        settings.setValue('exp_number', exp_num)

        # lattice parameters
        lattice_a = str(self.ui.lineEdit_a.text())
        settings.setValue('a', lattice_a)
        lattice_b = str(self.ui.lineEdit_b.text())
        settings.setValue('b', lattice_b)
        lattice_c = str(self.ui.lineEdit_c.text())
        settings.setValue('c', lattice_c)
        lattice_alpha = str(self.ui.lineEdit_alpha.text())
        settings.setValue('alpha', lattice_alpha)
        lattice_beta = str(self.ui.lineEdit_beta.text())
        settings.setValue('beta', lattice_beta)
        lattice_gamma = str(self.ui.lineEdit_gamma.text())
        settings.setValue('gamma', lattice_gamma)

        return

    def load_settings(self):
        """
        Load QSettings from previous saved file
        :return:
        """
        settings = QtCore.QSettings()

        # directories
        try:
            spice_dir = settings.value('local_spice_dir', '')
            self.ui.lineEdit_localSpiceDir.setText(str(spice_dir))
            work_dir = settings.value('work_dir')
            self.ui.lineEdit_workDir.setText(str(work_dir))

            # experiment number
            exp_num = settings.value('exp_number')
            self.ui.lineEdit_exp.setText(str(exp_num))

            # lattice parameters
            lattice_a = settings.value('a')
            self.ui.lineEdit_a.setText(str(lattice_a))
            lattice_b = settings.value('b')
            self.ui.lineEdit_b.setText(str(lattice_b))
            lattice_c = settings.value('c')
            self.ui.lineEdit_c.setText(str(lattice_c))
            lattice_alpha = settings.value('alpha')
            self.ui.lineEdit_alpha.setText(str(lattice_alpha))
            lattice_beta = settings.value('beta')
            self.ui.lineEdit_beta.setText(str(lattice_beta))
            lattice_gamma = settings.value('gamma')
            self.ui.lineEdit_gamma.setText(str(lattice_gamma))
        except TypeError as e:
            self.pop_one_button_dialog(str(e))
            return

        return

    def _get_lattice_parameters(self):
        """
        Get lattice parameters from GUI
        :return: (Boolean, Object).  True, 6-tuple as a, b, c, alpha, beta, gamm
                                     False: error message
        """
        status, ret_list = gutil.parse_float_editors([self.ui.lineEdit_a,
                                                      self.ui.lineEdit_b,
                                                      self.ui.lineEdit_c,
                                                      self.ui.lineEdit_alpha,
                                                      self.ui.lineEdit_beta,
                                                      self.ui.lineEdit_gamma])
        if status is False:
            err_msg = ret_list
            err_msg = 'Unable to parse unit cell due to %s' % err_msg
            return False, err_msg

        a, b, c, alpha, beta, gamma = ret_list

        return True, (a, b, c, alpha, beta, gamma)

    def _plot_raw_xml_2d(self, exp_no, scan_no, pt_no):
        """ Plot raw workspace from XML file for a measurement/pt.
        """
        # Check and load SPICE table file
        does_exist = self._myControl.does_spice_loaded(exp_no, scan_no)
        if does_exist is False:
            # Download data
            status, error_message = self._myControl.download_spice_file(exp_no, scan_no, over_write=False)
            if status is True:
                status, error_message = self._myControl.load_spice_scan_file(exp_no, scan_no)
                if status is False and self._allowDownload is False:
                    self.pop_one_button_dialog(error_message)
                    return
            else:
                self.pop_one_button_dialog(error_message)
                return
        # END-IF(does_exist)

        # Load Data for Pt's xml file
        does_exist = self._myControl.does_raw_loaded(exp_no, scan_no, pt_no)

        if does_exist is False:
            # Check whether needs to download
            status, error_message = self._myControl.download_spice_xml_file(scan_no, pt_no, exp_no=exp_no)
            if status is False:
                self.pop_one_button_dialog(error_message)
                return
            # Load SPICE xml file
            status, error_message = self._myControl.load_spice_xml_file(exp_no, scan_no, pt_no)
            if status is False:
                self.pop_one_button_dialog(error_message)
                return

        # Convert a list of vector to 2D numpy array for imshow()
        # Get data and plot
        raw_det_data = self._myControl.get_raw_detector_counts(exp_no, scan_no, pt_no)
        self.ui.graphicsView.clear_canvas()
        self.ui.graphicsView.add_plot_2d(raw_det_data, x_min=0, x_max=256, y_min=0, y_max=256,
                                         hold_prev_image=False)

        return
