#include "MantidQtMantidWidgets/DataProcessorUI/ParseNumerics.h"
#include <stdexcept>
#include "MantidQtMantidWidgets/WidgetDllOption.h"
namespace MantidQt {
namespace MantidWidgets {
double EXPORT_OPT_MANTIDQT_MANTIDWIDGETS parseDouble(QString const &in) {
  static auto ok = false;
  auto out = in.toDouble(&ok);
  if (ok)
    return out;
  else
    throw std::runtime_error("Failed to parse '" + in.toStdString() +
                             "' as a double numerical value.");
}

int EXPORT_OPT_MANTIDQT_MANTIDWIDGETS parseDenaryInteger(QString const &in) {
  static auto ok = false;
  auto constexpr BASE = 10;
  auto out = in.toInt(&ok, BASE);
  if (ok)
    return out;
  else
    throw std::runtime_error("Failed to parse '" + in.toStdString() +
                             "' as a denary integer.");
}
}
}
