
import unittest
import re
# Need to import mantid before we import SANSUtility
import mantid
from mantid.simpleapi import *
from mantid.api import mtd, WorkspaceGroup
from mantid.kernel import DateAndTime, time_duration, FloatTimeSeriesProperty,BoolTimeSeriesProperty,StringTimeSeriesProperty,StringPropertyWithValue
import SANSUtility as su
import re
import random

TEST_STRING_DATA = 'SANS2D0003434-add' + su.ADDED_EVENT_DATA_TAG
TEST_STRING_MON = 'SANS2D0003434-add_monitors' + su.ADDED_EVENT_DATA_TAG

TEST_STRING_DATA1 = TEST_STRING_DATA + '_1'
TEST_STRING_MON1 = TEST_STRING_MON + '_1'

TEST_STRING_DATA2 = TEST_STRING_DATA + '_2'
TEST_STRING_MON2 = TEST_STRING_MON + '_2'

TEST_STRING_DATA3 = TEST_STRING_DATA + '_3'
TEST_STRING_MON3 = TEST_STRING_MON + '_3'

def provide_group_workspace_for_added_event_data(event_ws_name, monitor_ws_name, out_ws_name):
    CreateWorkspace(DataX = [1,2,3], DataY = [2,3,4], OutputWorkspace = monitor_ws_name)
    CreateSampleWorkspace(WorkspaceType= 'Event', OutputWorkspace = event_ws_name)
    GroupWorkspaces(InputWorkspaces = [event_ws_name, monitor_ws_name ], OutputWorkspace = out_ws_name)

def addSampleLogEntry(log_name, ws, start_time, extra_time_shift):
    number_of_times = 10
    for i in range(0, number_of_times):

        val = random.randrange(0, 10, 1)
        date = DateAndTime(start_time)
        date +=  int(i*1e9)
        date += int(extra_time_shift*1e9)
        AddTimeSeriesLog(ws, Name=log_name, Time=date.__str__().strip(), Value=val)

def provide_event_ws_with_entries(name, start_time,number_events =0, extra_time_shift = 0.0, proton_charge = True, proton_charge_empty = False):
     # Create the event workspace
    ws = CreateSampleWorkspace(WorkspaceType= 'Event', NumEvents = number_events, OutputWorkspace = name)

    # Add the proton_charge log entries
    if proton_charge:
        if proton_charge_empty:
            addSampleLogEntry('proton_charge', ws, start_time, extra_time_shift)
        else:
            addSampleLogEntry('proton_charge', ws, start_time, extra_time_shift)

    # Add some other time series log entry
    addSampleLogEntry('time_series_2', ws, start_time, extra_time_shift)
    return ws

def provide_event_ws_custom(name, start_time, extra_time_shift = 0.0, proton_charge = True, proton_charge_empty = False):
    return provide_event_ws_with_entries(name=name, start_time=start_time,number_events = 100, extra_time_shift = extra_time_shift, proton_charge=proton_charge, proton_charge_empty=proton_charge_empty)

def provide_event_ws(name, start_time, extra_time_shift):
    return provide_event_ws_custom(name = name, start_time = start_time, extra_time_shift = extra_time_shift,  proton_charge = True)

def provide_event_ws_wo_proton_charge(name, start_time, extra_time_shift):
    return provide_event_ws_custom(name = name, start_time = start_time, extra_time_shift = extra_time_shift,  proton_charge = False)

# This test does not pass and was not used before 1/4/2015. SansUtilitytests was disabled.

#class TestSliceStringParser(unittest.TestCase):

#    def checkValues(self, list1, list2):

#        def _check_single_values( v1, v2):
#            self.assertAlmostEqual(v1, v2)

#        self.assertEqual(len(list1), len(list2))
#        for v1,v2 in zip(list1, list2):
#            start_1,stop_1 = v1
#            start_2, stop_2 = v2
#            _check_single_values(start_1, start_2)
#            _check_single_values(stop_1, stop_2)

#    def test_checkValues(self):
#        """sanity check to ensure that the others will work correctly"""
#        values = [  [[1,2],],
#                  [[None, 3],[4, None]],
#                 ]
#        for singlevalues in values:
#            self.checkValues(singlevalues, singlevalues)


#    def test_parse_strings(self):
#        inputs = { '1-2':[[1,2]],         # single period syntax  min < x < max
#                   '1.3-5.6':[[1.3,5.6]], # float
#                   '1-2,3-4':[[1,2],[3,4]],# more than one slice
#                   '>1':[[1, -1]],       # just lower bound
#                   '<5':[[-1, 5]],      # just upper bound
#                   '<5,8-9': [[-1, 5], [8,9]],
#                   '1:2:5': [[1,3], [3,5]] # sintax: start, step, stop
#            }

#        for (k, v) in inputs.items():
#            self.checkValues(su.sliceParser(k),v)

#    def test_accept_spaces(self):
#        self.checkValues(su.sliceParser("1 - 2, 3 - 4"), [[1,2],[3,4]])

#    def test_invalid_values_raise(self):
#        invalid_strs = ["5>6", ":3:", "MAX<min"]
#        for val in invalid_strs:
#            self.assertRaises(SyntaxError, su.sliceParser, val)

#    def test_empty_string_is_valid(self):
#        self.checkValues(su.sliceParser(""), [[-1,-1]])


class TestBundleAddedEventDataFilesToGroupWorkspaceFile(unittest.TestCase):
    def _prepare_workspaces(self, names):
        CreateSampleWorkspace(WorkspaceType = 'Event', OutputWorkspace = names[0])
        CreateWorkspace(DataX = [1,1,1], DataY = [1,1,1], OutputWorkspace = names[1])

        temp_save_dir = config['defaultsave.directory']
        if (temp_save_dir == ''):
            temp_save_dir = os.getcwd()

        data_file_name = os.path.join(temp_save_dir, names[0] + '.nxs')
        monitor_file_name = os.path.join(temp_save_dir, names[1] + '.nxs')

        SaveNexusProcessed(InputWorkspace = names[0], Filename = data_file_name)
        SaveNexusProcessed(InputWorkspace = names[1], Filename = monitor_file_name)

        file_names = [data_file_name, monitor_file_name]

        return file_names


    def _cleanup_workspaces(self, names):
        for name in names:
            DeleteWorkspace(name)


    def test_load_valid_added_event_data_and_monitor_file_produces_group_ws(self):
        # Arrange
        names = ['event_data', 'monitor']
        file_names = self._prepare_workspaces(names = names)
        self._cleanup_workspaces(names = names)

        # Act
        group_ws_name = 'g_ws'
        output_group_file_name = su.bundle_added_event_data_as_group(file_names[0], file_names[1])

        Load(Filename = output_group_file_name, OutputWorkspace = group_ws_name)
        group_ws = mtd[group_ws_name]

        # Assert
        self.assertTrue(isinstance(group_ws, WorkspaceGroup))
        self.assertEqual(group_ws.size(), 2)
        self.assertTrue(os.path.exists(file_names[0])) # File for group workspace exists
        self.assertFalse(os.path.exists(file_names[1]))  # File for monitors is deleted

        # Clean up
        ws_names_to_delete = []
        for ws_name in mtd.getObjectNames():
            if ws_name != group_ws_name:
                ws_names_to_delete.append(str(ws_name))
        self._cleanup_workspaces(names = ws_names_to_delete)

        if os.path.exists(file_names[0]):
            os.remove(file_names[0])

class TestLoadingAddedEventWorkspaceNameParsing(unittest.TestCase):

    def do_test_load_check(self, event_name, monitor_name):
        out_name = 'out_ws'
        provide_group_workspace_for_added_event_data(event_ws_name = event_name, monitor_ws_name = monitor_name, out_ws_name = out_name)
        out_ws = mtd[out_name]
        self.assertTrue(su.check_child_ws_for_name_and_type_for_added_eventdata(out_ws))
        DeleteWorkspace(out_ws)

    def test_check_regex_for_data(self):
        # Check when there is no special ending
        self.assertNotEqual(None, re.search(su.REG_DATA_NAME, TEST_STRING_DATA))
        # Check when there is a _1 ending
        self.assertNotEqual(None, re.search(su.REG_DATA_NAME, TEST_STRING_DATA1))
        # Check when there is a _2 ending
        self.assertNotEqual(None, re.search(su.REG_DATA_NAME, TEST_STRING_DATA2))
        # Check when there is a multiple ending
        self.assertNotEqual(None, re.search(su.REG_DATA_NAME, TEST_STRING_DATA3))


    def test_check_regex_for_data_monitors(self):
        # Check when there is no special ending
        self.assertNotEqual(None, re.search(su.REG_DATA_MONITORS_NAME, TEST_STRING_MON))
        # Check when there is a _1 ending
        self.assertNotEqual(None, re.search(su.REG_DATA_MONITORS_NAME, TEST_STRING_MON1))
        # Check when there is a _2 ending
        self.assertNotEqual(None, re.search(su.REG_DATA_MONITORS_NAME, TEST_STRING_MON2))
        # Check when there is a multiple ending
        self.assertNotEqual(None, re.search(su.REG_DATA_MONITORS_NAME, TEST_STRING_MON3))

    def test_regexes_do_not_clash(self):
        # Check when there is no special ending
        self.assertEqual(None, re.search(su.REG_DATA_NAME, TEST_STRING_MON)) 
        self.assertEqual(None, re.search(su.REG_DATA_MONITORS_NAME, TEST_STRING_DATA))
        # Check when there is a _1 ending
        self.assertEqual(None, re.search(su.REG_DATA_NAME, TEST_STRING_MON1)) 
        self.assertEqual(None, re.search(su.REG_DATA_MONITORS_NAME, TEST_STRING_DATA1))
        # Check when there is a _2 ending
        self.assertEqual(None, re.search(su.REG_DATA_NAME, TEST_STRING_MON2)) 
        self.assertEqual(None, re.search(su.REG_DATA_MONITORS_NAME, TEST_STRING_DATA2))
        # Check when there is a multiple ending
        self.assertEqual(None, re.search(su.REG_DATA_NAME, TEST_STRING_MON3)) 
        self.assertEqual(None, re.search(su.REG_DATA_MONITORS_NAME, TEST_STRING_DATA3))
    
    def test_check_child_file_names_for_valid_names(self):
        # Check when there is no special ending
        event_name = TEST_STRING_DATA
        monitor_name = TEST_STRING_MON
        self.do_test_load_check(event_name = event_name, monitor_name = monitor_name)

        # Check when there is a _1 ending
        event_name1 = TEST_STRING_DATA1
        monitor_name1 = TEST_STRING_MON1
        self.do_test_load_check(event_name = event_name1, monitor_name = monitor_name1)

        # Check when there is a _2 ending
        event_name2 = TEST_STRING_DATA2
        monitor_name2 = TEST_STRING_MON2
        self.do_test_load_check(event_name = event_name2, monitor_name = monitor_name2)

        # Check when there is a multiple ending
        event_name3 = TEST_STRING_DATA3
        monitor_name3 = TEST_STRING_MON3
        self.do_test_load_check(event_name = event_name3, monitor_name = monitor_name3)

class TestLoadingAddedEventWorkspaceExtraction(unittest.TestCase):
    _appendix = '_monitors'

    def do_test_extraction(self, event_name, monitor_name):
        out_ws_name = 'out_group'
        event_name_expect = out_ws_name
        monitor_name_expect = out_ws_name + self._appendix

        provide_group_workspace_for_added_event_data(event_ws_name = event_name, monitor_ws_name = monitor_name, out_ws_name = out_ws_name)
        out_ws_group = mtd[out_ws_name]
        
        # Act
        su.extract_child_ws_for_added_eventdata(out_ws_group, self._appendix)

        # Assert

        self.assertTrue(event_name_expect in mtd)
        self.assertTrue(monitor_name_expect in mtd)

        DeleteWorkspace(event_name_expect)
        DeleteWorkspace(monitor_name_expect)


    def test_extract_data_and_monitor_child_ws(self):
        # Check when there is no special ending
        self.do_test_extraction(TEST_STRING_DATA, TEST_STRING_MON)



class AddOperationTest(unittest.TestCase):
    def compare_added_workspaces(self,ws1, ws2, out_ws, start_time1, start_time2, extra_time_shift, isOverlay):
        self._compare_added_logs(ws1, ws2, out_ws, start_time1, start_time2, extra_time_shift, isOverlay)
        # Could compare events here, but trusting add algorithm for now

    def _compare_added_logs(self, ws1, ws2, out_ws,time1, time2, extra_time_shift, isOverlay):
        run1 = ws1.getRun()
        run2 = ws2.getRun()
        out_run = out_ws.getRun()
        props_out = out_run.getProperties()

        # Check that all times of workspace1 and workspace2 can be found in the outputworkspace
        for prop in props_out:
            if isinstance(prop, FloatTimeSeriesProperty) or isinstance(prop, BoolTimeSeriesProperty) or isinstance(prop, StringTimeSeriesProperty):
                prop1 = run1.getProperty(prop.name)
                prop2 = run2.getProperty(prop.name)
                self._compare_time_series(prop, prop1, prop2, time1, time2, extra_time_shift, isOverlay)
            elif isinstance(prop, StringPropertyWithValue):
                pass

    def _compare_time_series(self, prop_out, prop_in1, prop_in2, time1, time2, extra_time_shift, isOverlay):
        times_out = prop_out.times
        times1 = prop_in1.times
        times2 = prop_in2.times

        # Total time differnce is TIME1 - (TIME2 + extraShift)
        shift = 0.0
        if isOverlay:
            shift = time_duration.total_nanoseconds(DateAndTime(time1)- DateAndTime(time2))/1e9 - extra_time_shift

        # Check ws1 against output
        # We shift the second workspace onto the first workspace
        shift_ws1 = 0
        self._assert_times(times1, times_out, shift_ws1)
        # Check ws2 against output
        self._assert_times(times2, times_out, shift)
        # Check overlapping times
        self._compare_overlapping_times(prop_in1, prop_out, times1, times2, times_out, shift)

    def _assert_times(self, times_in, times_out, shift):
        for time in times_in:
            # Add the shift in nanaoseconds to the DateAndTime object ( the plus operator is defined
            # for nanoseconds)
            converted_time = time + int(shift*1e9)
            self.assertTrue(converted_time in times_out)

    def _compare_overlapping_times(self, prop_in1, prop_out, times1, times2, times_out, shift):
        overlap_times= []
        for time in times1:
            if time in times2:
                overlap_times.append(time)
        # Now go through all those overlap times and check that the value of the 
        # first workspace is recorded in the output
        for overlap_time in overlap_times:
            times1_list = list(times1)
            timesout_list = list(times_out)
            index_in1 = times1_list.index(overlap_time)
            index_out = timesout_list.index(overlap_time)
            value_in1 = prop_in1.nthValue(index_in1)
            value_out = prop_out.nthValue(index_out)
            self.assertTrue(value_in1 == value_out)

    def test_two_files_are_added_correctly_for_overlay_on(self):
        isOverlay = True
        names =['ws1', 'ws2']
        out_ws_name = 'out_ws'
        # Create event ws1
        start_time_1 = "2010-01-01T00:00:00"
        ws1 = provide_event_ws_with_entries(names[0],start_time_1, extra_time_shift = 0.0)
        # Create event ws2
        start_time_2 = "2012-01-01T00:10:00"
        ws2 = provide_event_ws(names[1],start_time_2, extra_time_shift = 0.0)
        # Create adder
        adder = su.AddOperation(isOverlay, '')
        # Act
        adder.add(ws1, ws2,out_ws_name, 0)
        out_ws = mtd[out_ws_name]
        # Assert
        self.compare_added_workspaces(ws1, ws2, out_ws, start_time_1, start_time_2, extra_time_shift = 0.0, isOverlay = isOverlay)

    def test_two_files_are_added_correctly_for_overlay_on_and_inverted_times(self):
        isOverlay = True
        names =['ws1', 'ws2']
        out_ws_name = 'out_ws'
        # Create event ws1
        start_time_1 = "2012-01-01T00:10:00"
        ws1 = provide_event_ws_with_entries(names[0],start_time_1, extra_time_shift = 0.0)
        # Create event ws2
        start_time_2 = "2010-01-01T00:00:00"
        ws2 = provide_event_ws(names[1],start_time_2, extra_time_shift = 0.0)
        # Create adder
        adder = su.AddOperation(isOverlay, '')
        # Act
        adder.add(ws1, ws2,out_ws_name, 0)
        out_ws = mtd[out_ws_name]
        # Assert
        self.compare_added_workspaces(ws1, ws2, out_ws, start_time_1, start_time_2, extra_time_shift = 0.0, isOverlay = isOverlay)

    def test_two_files_are_added_correctly_for_overlay_off(self):
        isOverlay = False
        names =['ws1', 'ws2']
        out_ws_name = 'out_ws'
        # Create event ws1
        start_time_1 = "2010-01-01T00:00:00"
        ws1 = provide_event_ws_with_entries(names[0],start_time_1, extra_time_shift = 0.0)
        # Create event ws2
        start_time_2 = "2012-01-01T01:00:00"
        ws2 = provide_event_ws(names[1],start_time_2, extra_time_shift = 0.0)
        # Create adder
        adder = su.AddOperation(False, '')
        # Act
        adder.add(ws1, ws2,out_ws_name, 0)
        out_ws = mtd[out_ws_name]
        # Assert
        self.compare_added_workspaces(ws1, ws2, out_ws, start_time_1, start_time_2, extra_time_shift = 0.0, isOverlay = isOverlay)

    def test_two_files_are_added_correctly_with_time_shift(self):
        isOverlay = True
        names =['ws1', 'ws2']
        out_ws_name = 'out_ws'
        time_shift = 100
        # Create event ws1
        start_time_1 = "2010-01-01T00:00:00"
        ws1 = provide_event_ws_with_entries(names[0],start_time_1, extra_time_shift = 0.0)
        # Create event ws2
        start_time_2 = "2012-01-01T01:10:00"
        ws2 = provide_event_ws(names[1],start_time_2, extra_time_shift = time_shift )
        # Create adder
        adder = su.AddOperation(True, '')
        # Act
        adder.add(ws1, ws2,out_ws_name, 0)
        out_ws = mtd[out_ws_name]
        # Assert
        self.compare_added_workspaces(ws1, ws2, out_ws, start_time_1, start_time_2, extra_time_shift = time_shift, isOverlay = isOverlay)


    def test_multiple_files_are_overlayed_correctly(self):
        isOverlay = True
        names =['ws1', 'ws2', 'ws3']
        out_ws_name = 'out_ws'
        out_ws_name2 = 'out_ws2'
        # Create event ws1
        start_time_1 = "2010-01-01T00:00:00"
        ws1 = provide_event_ws_with_entries(names[0],start_time_1, extra_time_shift = 0.0)
        # Create event ws2
        start_time_2 = "2012-01-01T00:00:00"
        ws2 = provide_event_ws(names[1],start_time_2, extra_time_shift = 0.0)
        # Create event ws3
        start_time_3 = "2013-01-01T00:00:00"
        ws3 = provide_event_ws(names[2],start_time_3, extra_time_shift = 0.0)
        # Create adder
        adder = su.AddOperation(True, '')
        # Act
        adder.add(ws1, ws2,out_ws_name, 0)
        adder.add(out_ws_name, ws3, out_ws_name2, 0)
        out_ws2 = mtd[out_ws_name2]
        out_ws = mtd[out_ws_name]
        # Assert
        self.compare_added_workspaces(out_ws, ws2, out_ws2, start_time_1, start_time_2, extra_time_shift = 0.0, isOverlay = isOverlay)

class TestCombineWorkspacesFactory(unittest.TestCase):
    def test_that_factory_returns_overlay_class(self):
        factory = su.CombineWorkspacesFactory()
        alg = factory.create_add_algorithm(True)
        self.assertTrue(isinstance(alg, su.OverlayWorkspaces))

    def test_that_factory_returns_overlay_class(self):
        factory = su.CombineWorkspacesFactory()
        alg = factory.create_add_algorithm(False)
        self.assertTrue(isinstance(alg, su.PlusWorkspaces))

class TestOverlayWorkspaces(unittest.TestCase):
    def test_time_from_proton_charge_log_is_recovered(self):
        # Arrange
        names =['ws1', 'ws2']
        out_ws_name = 'out_ws'

        start_time_1 = "2010-01-01T00:00:00"
        event_ws_1 = provide_event_ws(names[0],start_time_1, extra_time_shift = 0.0)

        start_time_2 = "2012-01-01T00:00:00"
        event_ws_2 = provide_event_ws(names[1],start_time_2, extra_time_shift = 0.0)

        # Act
        overlayWorkspaces = su.OverlayWorkspaces()
        time_difference = overlayWorkspaces._extract_time_difference_in_seconds(event_ws_1, event_ws_2)

        # Assert
        expected_time_difference = time_duration.total_nanoseconds(DateAndTime(start_time_1)- DateAndTime(start_time_2))/1e9
        self.assertEqual(time_difference, expected_time_difference)

        # Clean up 
        self._clean_up(names)
        self._clean_up(out_ws_name)

    def test_that_time_difference_adds_correct_optional_shift(self):
         # Arrange
        names =['ws1', 'ws2']
        out_ws_name = 'out_ws'

        start_time_1 = "2010-01-01T00:00:00"
        event_ws_1 = provide_event_ws(names[0],start_time_1, extra_time_shift = 0.0)

        # Extra shift in seconds
        optional_time_shift = 1000
        start_time_2 = "2012-01-01T00:00:00"
        event_ws_2 = provide_event_ws(names[1],start_time_2, extra_time_shift = optional_time_shift)

        # Act
        overlayWorkspaces = su.OverlayWorkspaces()
        time_difference = overlayWorkspaces._extract_time_difference_in_seconds(event_ws_1, event_ws_2)

        # Assert
        expected_time_difference = time_duration.total_nanoseconds(DateAndTime(start_time_1)- DateAndTime(start_time_2))/1e9
        expected_time_difference -= optional_time_shift # Need to subtract as we add the time shift to the subtrahend
        self.assertEqual(time_difference, expected_time_difference)

        # Clean up 
        self._clean_up(names)
        self._clean_up(out_ws_name)

    def test_error_is_raised_if_proton_charge_is_missing(self):
        # Arrange
        names =['ws1', 'ws2']
        out_ws_name = 'out_ws'

        start_time_1 = "2010-01-01T00:00:00"
        event_ws_1 = provide_event_ws_custom(name = names[0], start_time = start_time_1, extra_time_shift = 0.0, proton_charge = False)

        # Extra shift in seconds
        start_time_2 = "2012-01-01T00:00:00"
        event_ws_2 = provide_event_ws_custom(name = names[1], start_time = start_time_2, extra_time_shift = 0.0,proton_charge = False)

        # Act and Assert
        overlayWorkspaces = su.OverlayWorkspaces()
        with self.assertRaises(RuntimeError):
            overlayWorkspaces._extract_time_difference_in_seconds(event_ws_1, event_ws_2)

        # Clean up 
        self._clean_up(names)
        self._clean_up(out_ws_name)

    def test_correct_time_difference_is_extracted(self):
        pass
    def test_workspaces_are_normalized_by_proton_charge(self):
        pass

    def _clean_up(self,names):
        for name in names:
            if name in mtd.getObjectNames():
                DeleteWorkspace(name)

class TestTimeShifter(unittest.TestCase):
    def test_zero_shift_when_out_of_range(self):
        # Arrange
        time_shifts = ['12', '333.6', '-232']
        time_shifter = su.TimeShifter(time_shifts)

        # Act and Assert
        self.assertEqual(time_shifter.get_Nth_time_shift(0), 12.0)
        self.assertEqual(time_shifter.get_Nth_time_shift(1), 333.6)
        self.assertEqual(time_shifter.get_Nth_time_shift(2), -232.0)
        self.assertEqual(time_shifter.get_Nth_time_shift(3), 0.0)

    def test_zero_shift_when_bad_cast(self):
        # Arrange
        time_shifts = ['12', '33a.6', '-232']
        time_shifter = su.TimeShifter(time_shifts)

        # Act and Assert
        self.assertEqual(time_shifter.get_Nth_time_shift(0), 12.0)
        self.assertEqual(time_shifter.get_Nth_time_shift(1), 0.0)
        self.assertEqual(time_shifter.get_Nth_time_shift(2), -232.0)
        self.assertEqual(time_shifter.get_Nth_time_shift(3), 0.0)

if __name__ == "__main__":
    unittest.main()
