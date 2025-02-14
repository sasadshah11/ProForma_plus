#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/test_config.h>
#include <iostream>

using namespace OpenMS;

START_TEST(ProForma, "$Id$")

AASequence seq = AASequence::fromString("MGPEPTIDE");
ProForma pf(seq);

START_SECTION(registerCustomModification)
{
  pf.registerCustomModification("250.50", 250.50);
  pf.registerCustomModification("79.97", 79.97);
  pf.registerCustomModification("45.00", 45.00);

  TEST_EQUAL(ModificationsDB::getInstance()->has("250.50"), true);
  TEST_EQUAL(ModificationsDB::getInstance()->has("79.97"), true);
  TEST_EQUAL(ModificationsDB::getInstance()->has("45.00"), true);

  pf.registerCustomModification("250.50", 250.50);
  TEST_EQUAL(ModificationsDB::getInstance()->has("250.50"), true);
}
END_SECTION

START_SECTION(fromProFormaString)
{
  std::string proforma_str = "M[+250.50]PEPT[+79.97]ID[+45.00]E";
  AASequence result_seq = pf.fromProFormaString(proforma_str);

  TEST_STRING_EQUAL(result_seq.toString(), "M[+250.50]PEPT[+79.97]ID[+45.00]E");

  TEST_EXCEPTION(Exception::InvalidValue, pf.fromProFormaString("M[250.50]E"));
  TEST_EXCEPTION(Exception::InvalidValue, pf.fromProFormaString("M[+250.50E"));

  std::string invalid_brackets = "KHAESK(AD)[[Lys-loss,]]REEARAAVRGGK";
  TEST_EXCEPTION(Exception::InvalidValue, pf.fromProFormaString(invalid_brackets));

  std::string valid_brackets = "KHAESK(AD)[Lys-loss]REEARAAVRGGK";
  AASequence valid_seq = pf.fromProFormaString(valid_brackets);
  TEST_STRING_EQUAL(valid_seq.toString(), "KHAESK(AD)[Lys-loss]REEARAAVRGGK");
}
END_SECTION

START_SECTION(addModification_and_toProFormaString)
{
  pf.addModification(0, 0, "250.50", 250.50);
  pf.addModification(3, 3, "79.97", 79.97);
  pf.addModification(5, 5, "45.00", 45.00);

  String result = pf.toProFormaString();
  TEST_STRING_EQUAL(result, "M[+250.50]PEPT[+79.97]ID[+45.00]E");

  TEST_EXCEPTION(Exception::InvalidValue, pf.addModification(10, 10, "100.00", 100.00));
  TEST_EXCEPTION(Exception::InvalidValue, pf.addModification(0, 10, "100.00", 100.00));
  TEST_EXCEPTION(Exception::InvalidValue, pf.addModification(0, 0, "", 250.50));
}
END_SECTION

START_SECTION(handleInvalidModifications)
{
  TEST_EXCEPTION(Exception::InvalidValue, pf.addModification(0, 0, "InvalidModification", -999.99));
  TEST_EXCEPTION(Exception::InvalidValue, pf.addModification(0, 0, "InvalidModification", 1e10));
}
END_SECTION

START_SECTION(registerInvalidModifications)
{
  TEST_EXCEPTION(Exception::InvalidValue, pf.registerCustomModification("", 250.50));
  TEST_EXCEPTION(Exception::InvalidValue, pf.registerCustomModification("InvalidModification", -250.50));
}
END_SECTION

START_SECTION(boundaryConditions)
{
  ProForma pf_boundary(AASequence::fromString("A"));

  pf_boundary.addModification(0, 0, "100.00", 100.00);
  String result = pf_boundary.toProFormaString();
  TEST_STRING_EQUAL(result, "A[+100.00]");

  pf.addModification(0, 0, "50.00", 50.00);
  pf.addModification(8, 8, "30.00", 30.00);

  result = pf.toProFormaString();
  TEST_STRING_EQUAL(result, "M[+50.00]GPEPTID[+30.00]E");
}
END_SECTION

END_TEST
