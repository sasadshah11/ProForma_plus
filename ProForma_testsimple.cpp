#include <iostream>
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

void test_fixDoubleBracketsAndComma()
{
  OpenMS::AASequence seq = OpenMS::AASequence::fromString("KHAESKREEARAAVRGGK");
  OpenMS::ProForma proforma(seq);

  std::string proforma_str = "KHAESKREE(AR)[Lys-loss]AAVRGGK";

  proforma.fromProFormaString(proforma_str);
  std::string output = proforma.toProFormaString();
  if (output == proforma_str)
  {
    std::cout << "test_fixDoubleBracketsAndComma PASSED" << std::endl;
  }
  else
  {
    std::cerr << "test_fixDoubleBracketsAndComma FAILED: Expected " << proforma_str << " but got " << output << std::endl;
  }
}

void test_parseRangeModification()
{
  OpenMS::AASequence seq = OpenMS::AASequence::fromString("ACDEFGHIK");
  OpenMS::ProForma proforma(seq);

  std::string proforma_str = "A(CDE)[+12.50]FGHIK";
  proforma.fromProFormaString(proforma_str);
  std::string output = proforma.toProFormaString();
  if (output == proforma_str)
  {
    std::cout << "test_parseRangeModification PASSED" << std::endl;
  }
  else
  {
    std::cerr << "test_parseRangeModification FAILED: Expected " << proforma_str << " but got " << output << std::endl;
  }
}

void test_parseSimpleModification()
{
  OpenMS::AASequence seq = OpenMS::AASequence::fromString("ACDEFGHIK");
  OpenMS::ProForma proforma(seq);

  std::string proforma_str = "A[Phospho]CDEFGHIK";
  proforma.fromProFormaString(proforma_str);
  std::string output = proforma.toProFormaString();
  if (output == proforma_str)
  {
    std::cout << "test_parseSimpleModification PASSED" << std::endl;
  }
  else
  {
    std::cerr << "test_parseSimpleModification FAILED: Expected " << proforma_str << " but got " << output << std::endl;
  }
}
void test_parseMassShiftModification()
{
  OpenMS::AASequence seq = OpenMS::AASequence::fromString("ACDEFGHIK");
  OpenMS::ProForma proforma(seq);

  std::string proforma_str = "A[+15.99]CDEFGHIK";
  proforma.fromProFormaString(proforma_str);
  std::string output = proforma.toProFormaString();
  if (output == proforma_str)
  {
    std::cout << "test_parseMassShiftModification PASSED" << std::endl;
  }
  else
  {
    std::cerr << "test_parseMassShiftModification FAILED: Expected " << proforma_str << " but got " << output << std::endl;
  }
}

void test_parseNTerminalModification()
{
  OpenMS::AASequence seq = OpenMS::AASequence::fromString("ACDEFGHIK");
  OpenMS::ProForma proforma(seq);

  std::string proforma_str = "[Acetyl]-ACDEFGHIK";
  proforma.fromProFormaString(proforma_str);
  std::string output = proforma.toProFormaString();
  if (output == proforma_str)
  {
    std::cout << "test_parseNTerminalModification PASSED" << std::endl;
  }
  else
  {
    std::cerr << "test_parseNTerminalModification FAILED: Expected " << proforma_str << " but got " << output << std::endl;
  }
}

void test_parseCTerminalModification()
{
  OpenMS::AASequence seq = OpenMS::AASequence::fromString("ACDEFGHIK");
  OpenMS::ProForma proforma(seq);

  std::string proforma_str = "ACDEFGHIK-[Amidation]";
  proforma.fromProFormaString(proforma_str);
  std::string output = proforma.toProFormaString();
  if (output == proforma_str)
  {
    std::cout << "test_parseCTerminalModification PASSED" << std::endl;
  }
  else
  {
    std::cerr << "test_parseCTerminalModification FAILED: Expected " << proforma_str << " but got " << output << std::endl;
  }
}

int main()
{
  test_fixDoubleBracketsAndComma();
  test_parseRangeModification();
  test_parseSimpleModification();
  test_parseMassShiftModification();
  test_parseNTerminalModification();
  test_parseCTerminalModification();
  return 0;
}
