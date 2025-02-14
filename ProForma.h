// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Ayesha Feroz $
// $Authors: Ayesha Feroz, Tom Müller$
// --------------------------------------------------------------------------
#ifndef OPENMS_CHEMISTRY_PROFORMA_H
#define OPENMS_CHEMISTRY_PROFORMA_H

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace OpenMS
{
struct ModificationAttributes
{
  double mass_shift = 0.0;
  std::string modification_name;
  size_t end_pos = 0;
};

class OPENMS_DLLAPI ProForma
{
public:
  explicit ProForma(const AASequence& seq);

  AASequence fromProFormaString(const std::string& proforma_str);
  std::string toProFormaString() const;

  void addModification(size_t start_pos, size_t end_pos, const std::string& mod_id, double mass_shift);
  void registerCustomModification(const String& id, double mass_shift);

private:
  AASequence sequence_;
  std::unordered_map<size_t, ModificationAttributes> modifications_;
  std::unordered_set<std::string> supported_cvs_{"UNIMOD", "MOD", "RESID", "XLMOD", "GNO"};

  void parseModifications(const std::string& modString, size_t& pos, size_t residue_pos);
  void throwParseError(const std::string& message) const;
  void validateCVModification(const std::string& modification);
};

} // namespace OpenMS

#endif // OPENMS_CHEMISTRY_PROFORMA_H
