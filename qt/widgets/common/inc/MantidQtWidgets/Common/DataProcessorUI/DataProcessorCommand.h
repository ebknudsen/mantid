#ifndef MANTIDQTMANTIDWIDGETS_DATAPROCESSORCOMMAND_H
#define MANTIDQTMANTIDWIDGETS_DATAPROCESSORCOMMAND_H

#include <memory>
#include <vector>

#include <QString>

namespace MantidQt {
namespace MantidWidgets {
/** @class DataProcessorCommand

DataProcessorCommand is an interface which defines the functions any data
processor action needs to support.

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
class DataProcessorCommand {
public:
  DataProcessorCommand() : m_child(){};
  DataProcessorCommand(DataProcessorCommand&&) = delete;
  DataProcessorCommand &operator=(DataProcessorCommand&&) = delete;
  virtual ~DataProcessorCommand(){};


  virtual void execute() = 0;
  virtual QString name() const = 0;
  virtual QString icon() const = 0;
  virtual QString tooltip() const = 0;
  virtual QString whatsthis() const = 0;
  virtual QString shortcut() const = 0;
  virtual bool modifiesTable() const = 0;
  virtual bool hasChild() const final { return !m_child.empty(); };
  virtual void
  setChild(std::vector<std::unique_ptr<DataProcessorCommand>> child) final {
    m_child = std::move(child);
  }
  virtual std::vector<std::unique_ptr<DataProcessorCommand>> &getChild() final {
    return m_child;
  }
  virtual bool isSeparator() const final {
    return name().isEmpty() && icon().isEmpty();
  }

protected:
  std::vector<std::unique_ptr<DataProcessorCommand>> m_child;
};

/// Typedef for a shared pointer to \c ReflSearchModel
typedef std::unique_ptr<DataProcessorCommand> DataProcessorCommand_uptr;
}
}
#endif /*MANTIDQTMANTIDWIDGETS_DATAPROCESSORCOMMAND_H*/
