#include "MantidGeometry/Crystal/PointGroup.h"
#include "MantidKernel/System.h"

#include <set>
#include <boost/make_shared.hpp>
#include <iostream>

#include "MantidGeometry/Crystal/PointGroupFactory.h"
#include "MantidGeometry/Crystal/SymmetryOperationFactory.h"

namespace Mantid {
namespace Geometry {
using Kernel::V3D;
using Kernel::IntMatrix;

/**
 * Returns all equivalent reflections for the supplied hkl.
 *
 * This method returns a vector containing all equivalent hkls for the supplied
 * one. It depends on the internal state of the pointgroup object (e.g. which
 * symmetry operations and therefore, which transformation matrices are
 * present). This internal state is unique for each concrete point group and
 * is set in the constructor.
 *
 * The returned vector always contains a set of unique hkls, so for special hkls
 * like (100), it has fewer entries than for a general hkl. See also
 * PointGroup::getEquivalentSet.
 *
 * @param hkl :: Arbitrary hkl
 * @return :: std::vector containing all equivalent hkls.
 */
std::vector<V3D> PointGroup::getEquivalents(const V3D &hkl) const {
  return getEquivalentSet(hkl);

  // return std::vector<V3D>(equivalents.rbegin(), equivalents.rend());
}

/**
 * Returns the same V3D for all equivalent hkls.
 *
 * This method is closely related to PointGroup::getEquivalents. It returns the
 * same V3D for all hkls of one "family". For example in a cubic point group
 * it will return (100) for (001), (010), (0-10), etc.
 *
 * It can be used to generate a set of symmetry independent hkls, useful for
 * example in powder diffraction.
 *
 * @param hkl :: Arbitrary hkl
 * @return :: hkl specific to a family of index-triplets
 */
V3D PointGroup::getReflectionFamily(const Kernel::V3D &hkl) const {
  return *getEquivalentSet(hkl).rbegin();
}

bool PointGroup::groupHasNoTranslations(const Group &group) const {
  const std::vector<SymmetryOperation> &symOps = group.getSymmetryOperations();

  for (auto op = symOps.begin(); op != symOps.end(); ++op) {
    if ((*op).hasTranslation()) {
      return false;
    }
  }

  return true;
}

/// Protected constructor - can not be used directly.
PointGroup::PointGroup(const std::string &symbolHM, const Group &group,
                       const std::string &name)
    : Group(group), m_permutations(), m_symbolHM(symbolHM), m_name(name) {
  createPermutationsFromOperations();
}

PointGroup::PointGroup(const PointGroup &other)
    : Group(other), m_permutations(), m_symbolHM(other.m_symbolHM),
      m_name(other.m_name) {
  createPermutationsFromOperations();
}

PointGroup &PointGroup::operator=(const PointGroup &other) {
  Group::operator=(other);

  createPermutationsFromOperations();

  m_symbolHM = other.m_symbolHM;
  m_name = other.m_name;

  return *this;
}

/// Hermann-Mauguin symbol
std::string PointGroup::getSymbol() const { return m_symbolHM; }

bool PointGroup::isEquivalent(const Kernel::V3D &hkl,
                              const Kernel::V3D &hkl2) const {
  std::vector<V3D> hklEquivalents = getEquivalentSet(hkl);

  return (std::find(hklEquivalents.begin(), hklEquivalents.end(), hkl2) !=
          hklEquivalents.end());
}

/**
 * Generates a set of hkls
 *
 * This method applies all transformation matrices to the supplied hkl and puts
 * it into a set, which is returned in the end. Using a set ensures that each
 * hkl occurs once and only once. This set is the set of equivalent hkls,
 * specific to a concrete point group.
 *
 * The symmetry operations need to be set prior to calling this method by a call
 * to PointGroup::setTransformationMatrices.
 *
 * @param hkl :: Arbitrary hkl
 * @return :: set of hkls.
 */
std::vector<V3D> PointGroup::getEquivalentSet(const Kernel::V3D &hkl) const {
  std::vector<V3D> equivalents;
  equivalents.reserve(m_permutations.size());

  for (auto op = m_permutations.begin(); op != m_permutations.end(); ++op) {
    equivalents.push_back((*op).getPermutation(hkl));
  }

  std::sort(equivalents.begin(), equivalents.end());

  equivalents.erase(std::unique(equivalents.begin(), equivalents.end()),
                    equivalents.end());

  return equivalents;
}

void PointGroup::createPermutationsFromOperations() {
  m_permutations.clear();

  const std::vector<SymmetryOperation> &symOps = getSymmetryOperations();
  for (auto op = symOps.begin(); op != symOps.end(); ++op) {
    m_permutations.push_back((*op).matrix());
  }
}

/** @return a vector with all possible PointGroup objects */
std::vector<PointGroup_sptr> getAllPointGroups() {
  std::vector<std::string> allSymbols =
      PointGroupFactory::Instance().getAllPointGroupSymbols();

  std::vector<PointGroup_sptr> out;
  for (auto it = allSymbols.begin(); it != allSymbols.end(); ++it) {
    out.push_back(PointGroupFactory::Instance().createPointGroup(*it));
  }

  return out;
}

PointGroupCrystalSystemMap getPointGroupsByCrystalSystem() {
  PointGroupCrystalSystemMap map;

  std::vector<PointGroup_sptr> pointGroups = getAllPointGroups();
  for (size_t i = 0; i < pointGroups.size(); ++i) {
    map.insert(std::make_pair(pointGroups[i]->crystalSystem(), pointGroups[i]));
  }

  return map;
}

} // namespace Mantid
} // namespace Geometry
