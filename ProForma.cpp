// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Ayesha Feroz $
// $Authors: Ayesha Feroz, Tom Müller$
// --------------------------------------------------------------------------
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <iostream>

namespace OpenMS
{

ProForma::ProForma(const AASequence& seq) : sequence_(seq) {}

void ProForma::throwParseError(const std::string& message) const
{
  throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message, "Invalid ProForma input");
}

void ProForma::validateCVModification(const std::string& modification)
{
  size_t colon_pos = modification.find(':');
  if (colon_pos == std::string::npos)
    throwParseError("No CV prefix found in modification: " + modification);

  std::string cv = modification.substr(0, colon_pos);
  if (supported_cvs_.find(cv) == supported_cvs_.end())
    throwParseError("Unsupported CV/ontology: " + cv);
}

void ProForma::registerCustomModification(const String& id, double mass_shift)
{
  if (id.empty() || mass_shift == 0)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Invalid modification ID or mass shift: " + static_cast<std::string>(id), "InvalidModification");
  }

  if (!ModificationsDB::getInstance()->has(id))
  {
    auto custom_mod = std::make_unique<ResidueModification>();
    custom_mod->setId(id);
    custom_mod->setName("Custom Modification");
    custom_mod->setFullId("Custom|" + id);
    custom_mod->setDiffMonoMass(mass_shift);
    custom_mod->setTermSpecificity(ResidueModification::ANYWHERE);
    ModificationsDB::getInstance()->addModification(std::move(custom_mod));
  }
}

void ProForma::parseModifications(const std::string& modString, size_t& pos, size_t residue_pos)
{
  size_t modStart = modString.find('[', pos);
  size_t modEnd = modString.find(']', modStart);
  if (modStart == std::string::npos || modEnd == std::string::npos)
    throwParseError("Invalid modification format.");

  std::string mod_str = modString.substr(modStart + 1, modEnd - modStart - 1);

  if (mod_str.empty() || (mod_str[0] != '+' && mod_str[0] != '-' && mod_str.find("Cation:") == std::string::npos))
  {
    throwParseError("Invalid ProForma input. Modification must start with '+' or '-' or be a recognized cation.");
  }

  try
  {
    double mass_shift = std::stod(mod_str);
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(2) << mass_shift;
    std::string custom_mod_id = ss.str();

    registerCustomModification(custom_mod_id, mass_shift);
    modifications_[residue_pos - 1] = {mass_shift, custom_mod_id, residue_pos - 1};
  }
  catch (const std::exception& e)
  {
    throwParseError("Error parsing modification: " + mod_str);
  }

  pos = modEnd + 1;
}

AASequence ProForma::fromProFormaString(const std::string& proforma_str)
{
  std::string sequence_str;
  size_t pos = 0;
  size_t residue_pos = 0;

  while (pos < proforma_str.size())
  {
    if (std::isalpha(proforma_str[pos]))
    {
      sequence_str += proforma_str[pos];
      residue_pos++;
      pos++;
    }
    else if (proforma_str[pos] == '[')
    {
      parseModifications(proforma_str, pos, residue_pos);
    }
    else
    {
      pos++;
    }
  }

  AASequence seq = AASequence::fromString(sequence_str);

  for (const auto& mod : modifications_)
  {
    size_t position = mod.first;
    const ModificationAttributes& attributes = mod.second;

    if (position < seq.size() && !attributes.modification_name.empty())
    {
      seq.setModification(position, attributes.modification_name);
    }
    else
    {
      std::cerr << "Error: Invalid modification at position " << position << std::endl;
    }
  }

  return seq;
}

void ProForma::addModification(size_t start_pos, size_t end_pos, const std::string& mod_id, double mass_shift)
{
  if (start_pos >= sequence_.size() || end_pos >= sequence_.size() || start_pos > end_pos)
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Invalid start or end positions for modification.", "Invalid Position");
  }

  for (size_t pos = start_pos; pos <= end_pos; ++pos)
  {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(2) << mass_shift;
    std::string custom_mod_id = mod_id.empty() ? ss.str() : mod_id;

    registerCustomModification(custom_mod_id, mass_shift);
    modifications_[pos] = {mass_shift, custom_mod_id, end_pos};
  }
}

std::string ProForma::toProFormaString() const
{
  std::ostringstream ss;
  for (size_t i = 0; i < sequence_.size(); ++i)
  {
    ss << sequence_.getResidue(i).getOneLetterCode();
    auto it = modifications_.find(i);

    if (it != modifications_.end())
    {
      const ModificationAttributes& attr = it->second;
      if (!attr.modification_name.empty())
      {
        ss << "[" << (attr.mass_shift > 0 ? "+" : "")
           << std::fixed << std::setprecision(2) << attr.mass_shift << "]";
      }
    }
  }
  return ss.str();
}

}
