#ifndef MANTID_MDEVENTS_COORDTRANSFORMPARSERTEST_H_
#define MANTID_MDEVENTS_COORDTRANSFORMPARSERTEST_H_

#include "MantidMDEvents/CoordTransformParser.h"
#include "MantidMDEvents/CoordTransform.h"
#include <cxxtest/TestSuite.h>

#include <Poco/DOM/DOMParser.h>
#include <Poco/DOM/Document.h>
#include <Poco/DOM/Element.h>
#include <Poco/DOM/NodeList.h>
#include <Poco/DOM/NodeIterator.h>
#include <Poco/DOM/NodeFilter.h>
#include <Poco/File.h>
#include <Poco/Path.h>

using namespace Mantid::MDEvents;

class CoordTransformParserTest : public CxxTest::TestSuite
{
private:

  class MockCoordTransformParser : public CoordTransformParser
  {
    virtual CoordTransform* createTransform(Poco::XML::Element* coordTransElement) const
    {
      return new CoordTransform(1, 1);
    }
  };

public:

  void testSuccessfulParse()
  {
    std::string xmlToParse = std::string("<CoordTransform>") +
    "<Type>CoordTransform</Type>" +
    "<ParameterList>" +
    "<Parameter><Type>InDimParameter</Type><Value>2</Value></Parameter>" + 
    "<Parameter><Type>OutDimParameter</Type><Value>2</Value></Parameter>" + 
    "<Parameter><Type>AffineMatrixParameter</Type><Value>0,1,2;3,4,5;6,7,8</Value></Parameter>" + 
    "</ParameterList></CoordTransform>";

   Poco::XML::DOMParser pParser;
   Poco::XML::Document* pDoc = pParser.parseString(xmlToParse);
   Poco::XML::Element* pRootElem = pDoc->documentElement();

   CoordTransformParser parser;
   CoordTransform* transform = parser.createTransform(pRootElem);

   AffineMatrixType product = transform->getMatrix();

   //Check that matrix is recovered.
   TS_ASSERT_EQUALS(0, product[0][0]);
   TS_ASSERT_EQUALS(1, product[0][1]);
   TS_ASSERT_EQUALS(2, product[0][2]);
   TS_ASSERT_EQUALS(3, product[1][0]);
   TS_ASSERT_EQUALS(4, product[1][1]);
   TS_ASSERT_EQUALS(5, product[1][2]);
   TS_ASSERT_EQUALS(6, product[2][0]);
   TS_ASSERT_EQUALS(7, product[2][1]);
   TS_ASSERT_EQUALS(8, product[2][2]);

   //Circular check. Acutally hard to debug, but gives certainty that serialization and deserialization cause no side effects.
   TSM_ASSERT_EQUALS("Parsing has not occured correctly if the output is not equal to the intput", transform->toXMLString(), xmlToParse);

   delete transform;

  }

  void testNotACoordTransformThrows()
  {
   std::string xmlToParse = std::string("<OTHER></OTHER>");

   Poco::XML::DOMParser pParser;
   Poco::XML::Document* pDoc = pParser.parseString(xmlToParse);
   Poco::XML::Element* pRootElem = pDoc->documentElement();

   CoordTransformParser parser;
   TSM_ASSERT_THROWS("XML root node must be a coordinate transform", parser.createTransform(pRootElem), std::invalid_argument);
  }

  void testNoSuccessorThrows()
  {
   std::string xmlToParse = "<CoordTransform><Type>OTHER</Type></CoordTransform>"; //type is not a coordinate transform, so should try to use it's successor

   Poco::XML::DOMParser pParser;
   Poco::XML::Document* pDoc = pParser.parseString(xmlToParse);
   Poco::XML::Element* pRootElem = pDoc->documentElement();

   CoordTransformParser parser;
   TSM_ASSERT_THROWS("Should throw since no successor parser has been set", parser.createTransform(pRootElem), std::runtime_error);
  }

  void testDelegateToSuccessor()
  {
   std::string xmlToParse = "<CoordTransform><Type>OTHER</Type></CoordTransform>"; //type is not a coordinate transform, so should try to use it's successor

   Poco::XML::DOMParser pParser;
   Poco::XML::Document* pDoc = pParser.parseString(xmlToParse);
   Poco::XML::Element* pRootElem = pDoc->documentElement();

   CoordTransformParser parser;
   parser.setSuccessor(new MockCoordTransformParser);
   CoordTransform* product = parser.createTransform(pRootElem);
   delete product;
  }

};

#endif