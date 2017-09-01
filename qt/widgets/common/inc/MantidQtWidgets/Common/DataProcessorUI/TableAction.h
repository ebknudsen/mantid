#ifndef MANTID_CUSTOMINTERFACES_REFLACTION_H
#define MANTID_CUSTOMINTERFACES_REFLACTION_H
namespace MantidQt {
namespace MantidWidgets {
/** @enum TableAction

TableAction is an enumeration of actions concerned with the entire table.

Copyright &copy; 2011-16 ISIS Rutherford Appleton Laboratory, NScD Oak Ridge
National Laboratory & European Spallation Source

This file is part of Mantid.

Mantid is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

Mantid is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

File change history is stored at: <https://github.com/mantidproject/mantid>.
Code Documentation is available at: <http://doxygen.mantidproject.org>
*/
enum class TableAction {
  OPEN_TABLE,
  NEW_TABLE,
  SAVE_TABLE,
  SAVE_TABLE_AS,
  IMPORT_TBL_FILE,
  EXPORT_TBL_FILE,
  OPTIONS
};
}
}
#endif /*MANTID_CUSTOMINTERFACES_REFLACTION_H*/