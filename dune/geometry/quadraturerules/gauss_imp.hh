// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
//
// WARNING
// This file is automatically generated by jacobian.mac! Don't edit by hand!

#ifndef DUNE_INCLUDING_IMPLEMENTATION
#error This is a private header that should not be included directly.
#error Use #include <dune/geometry/quadraturerules.hh> instead.
#endif
#undef DUNE_INCLUDING_IMPLEMENTATION

namespace Dune {

  // for fundamental types
  template<typename ct>
  void GaussQuadratureInitHelper<ct,true>::init(int p,
         std::vector< FieldVector<ct, 1> > & _points,
         std::vector< ct > & _weight,
         int & delivered_order)
  {
    switch(p)
    {
    // order 0,1
    case 0 :
    case 1 :
      delivered_order = 1;
      _points.resize(1);
      _weight.resize(1);
      _points[0] = 0.5;
      _weight[0] = 1.0;
      break;

    // order 2,3
    case 2 :
    case 3 :
      delivered_order = 3;
      _points.resize(2);
      _weight.resize(2);
      _points[0] = 0.2113248654051871177454256097490212721761991243649365619906988367580111638485333271531423022071252374;
      _weight[0] = 0.5;
      _points[1] = 0.7886751345948128822545743902509787278238008756350634380093011632419888361514666728468576977928747626;
      _weight[1] = 0.5;
      break;

    // order 4,5
    case 4 :
    case 5 :
      delivered_order = 5;
      _points.resize(3);
      _weight.resize(3);
      _points[0] = 0.8872983346207416885179265399782399610832921705291590826587573766113483091936979033519287376858673518;
      _weight[0] = 0.2777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777;
      _points[1] = 0.1127016653792583114820734600217600389167078294708409173412426233886516908063020966480712623141326482;
      _weight[1] = 0.2777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777778;
      _points[2] = 0.5;
      _weight[2] = 0.4444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444;
      break;

    // order 6,7
    case 6 :
    case 7 :
      delivered_order = 7;
      _points.resize(4);
      _weight.resize(4);
      _points[0] = 0.9305681557970262876119732444464047525478626898148588188078609604532647357475244328520811699422396526;
      _weight[0] = 0.1739274225687269286865319746109997036176743479169467702462646597593759337329551758609918838661290797;
      _points[1] = 0.0694318442029737123880267555535952474521373101851411811921390395467352642524755671479188300577603474;
      _weight[1] = 0.1739274225687269286865319746109997036176743479169467702462646597593759337329551758609918838661290798;
      _points[2] = 0.3300094782075718675986671204483776563997120651145428237035230115894899847683814827610623597822225942;
      _weight[2] = 0.3260725774312730713134680253890002963823256520830532297537353402406240662670448241390081161338709202;
      _points[3] = 0.6699905217924281324013328795516223436002879348854571762964769884105100152316185172389376402177774058;
      _weight[3] = 0.3260725774312730713134680253890002963823256520830532297537353402406240662670448241390081161338709202;
      break;

    // order 8,9
    case 8 :
    case 9 :
      delivered_order = 9;
      _points.resize(5);
      _weight.resize(5);
      _points[0] = 0.953089922969331996398813439149696482562825955381265431436881143271885397458343423470571494776771131;
      _weight[0] = 0.1184634425280945437571320203599586813216300011062070077914139441108586442015215492899967152469757223;
      _points[1] = 0.04691007703066800360118656085030351743717404461873456856311885672811460254165657652942850522322886904;
      _weight[1] = 0.1184634425280945437571320203599586813216300011062070077914139441108586442015215492899967152469757223;
      _points[2] = 0.7692346550528415455181572103501044024836433034527799781011158135297355926838776455179018336252854658;
      _weight[2] = 0.2393143352496832340206457574178190964561477766715707699863638336669191335762562284877810625308020554;
      _points[3] = 0.2307653449471584544818427896498955975163566965472200218988841864702644073161223544820981663747145342;
      _weight[3] = 0.2393143352496832340206457574178190964561477766715707699863638336669191335762562284877810625308020554;
      _points[4] = 0.5;
      _weight[4] = 0.2844444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444;
      break;

    // order 10,11
    case 10 :
    case 11 :
      delivered_order = 11;
      _points.resize(6);
      _weight.resize(6);
      _points[0] = 0.03376524289842398609384922275300269543261713114385508756372519173669324957789990186185563003903700748;
      _weight[0] = 0.08566224618958517252014807108636644676341125074202199119931771989947288027117007732396385271319433444;
      _points[1] = 0.9662347571015760139061507772469973045673828688561449124362748082633067504221000981381443699609629925;
      _weight[1] = 0.08566224618958517252014807108636644676341125074202199119931771989947288027117007732396385271319433508;
      _points[2] = 0.8306046932331322568306997975099526735032242821975850354072633529260917483035715504721432018732307282;
      _weight[2] = 0.1803807865240693037849167569188580558307609463733727411448696201185700189186308591604811009944096738;
      _points[3] = 0.1693953067668677431693002024900473264967757178024149645927366470739082516964284495278567981267692718;
      _weight[3] = 0.180380786524069303784916756918858055830760946373372741144869620118570018918630859160481100994409674;
      _points[4] = 0.6193095930415984543152508608403559677093053150700106750906975822871374671378199211246122136286745658;
      _weight[4] = 0.2339569672863455236949351719947754974058278028846052676558126599819571008101990635155550462923959912;
      _points[5] = 0.3806904069584015456847491391596440322906946849299893249093024177128625328621800788753877863713254342;
      _weight[5] = 0.2339569672863455236949351719947754974058278028846052676558126599819571008101990635155550462923959912;
      break;

    default :
      DUNE_THROW(QuadratureOrderOutOfRange, "Quadrature rule " << p << " not supported!");
    }
  }

  // for non-fundamental types: assign numbers as strings
  template<typename ct>
  void GaussQuadratureInitHelper<ct,false>::init(int p,
         std::vector< FieldVector<ct, 1> > & _points,
         std::vector< ct > & _weight,
         int & delivered_order)
  {
    switch(p)
    {
    // order 0,1
    case 0 :
    case 1 :
      delivered_order = 1;
      _points.resize(1);
      _weight.resize(1);
      _points[0] = ct("0.5");
      _weight[0] = ct("1.0");
      break;

    // order 2,3
    case 2 :
    case 3 :
      delivered_order = 3;
      _points.resize(2);
      _weight.resize(2);
      _points[0] = ct("0.2113248654051871177454256097490212721761991243649365619906988367580111638485333271531423022071252374");
      _weight[0] = ct("0.5");
      _points[1] = ct("0.7886751345948128822545743902509787278238008756350634380093011632419888361514666728468576977928747626");
      _weight[1] = ct("0.5");
      break;

    // order 4,5
    case 4 :
    case 5 :
      delivered_order = 5;
      _points.resize(3);
      _weight.resize(3);
      _points[0] = ct("0.8872983346207416885179265399782399610832921705291590826587573766113483091936979033519287376858673518");
      _weight[0] = ct("0.2777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777");
      _points[1] = ct("0.1127016653792583114820734600217600389167078294708409173412426233886516908063020966480712623141326482");
      _weight[1] = ct("0.2777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777778");
      _points[2] = ct("0.5");
      _weight[2] = ct("0.4444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444");
      break;

    // order 6,7
    case 6 :
    case 7 :
      delivered_order = 7;
      _points.resize(4);
      _weight.resize(4);
      _points[0] = ct("0.9305681557970262876119732444464047525478626898148588188078609604532647357475244328520811699422396526");
      _weight[0] = ct("0.1739274225687269286865319746109997036176743479169467702462646597593759337329551758609918838661290797");
      _points[1] = ct("0.0694318442029737123880267555535952474521373101851411811921390395467352642524755671479188300577603474");
      _weight[1] = ct("0.1739274225687269286865319746109997036176743479169467702462646597593759337329551758609918838661290798");
      _points[2] = ct("0.3300094782075718675986671204483776563997120651145428237035230115894899847683814827610623597822225942");
      _weight[2] = ct("0.3260725774312730713134680253890002963823256520830532297537353402406240662670448241390081161338709202");
      _points[3] = ct("0.6699905217924281324013328795516223436002879348854571762964769884105100152316185172389376402177774058");
      _weight[3] = ct("0.3260725774312730713134680253890002963823256520830532297537353402406240662670448241390081161338709202");
      break;

    // order 8,9
    case 8 :
    case 9 :
      delivered_order = 9;
      _points.resize(5);
      _weight.resize(5);
      _points[0] = ct("0.953089922969331996398813439149696482562825955381265431436881143271885397458343423470571494776771131");
      _weight[0] = ct("0.1184634425280945437571320203599586813216300011062070077914139441108586442015215492899967152469757223");
      _points[1] = ct("0.04691007703066800360118656085030351743717404461873456856311885672811460254165657652942850522322886904");
      _weight[1] = ct("0.1184634425280945437571320203599586813216300011062070077914139441108586442015215492899967152469757223");
      _points[2] = ct("0.7692346550528415455181572103501044024836433034527799781011158135297355926838776455179018336252854658");
      _weight[2] = ct("0.2393143352496832340206457574178190964561477766715707699863638336669191335762562284877810625308020554");
      _points[3] = ct("0.2307653449471584544818427896498955975163566965472200218988841864702644073161223544820981663747145342");
      _weight[3] = ct("0.2393143352496832340206457574178190964561477766715707699863638336669191335762562284877810625308020554");
      _points[4] = ct("0.5");
      _weight[4] = ct("0.2844444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444");
      break;

    // order 10,11
    case 10 :
    case 11 :
      delivered_order = 11;
      _points.resize(6);
      _weight.resize(6);
      _points[0] = ct("0.03376524289842398609384922275300269543261713114385508756372519173669324957789990186185563003903700748");
      _weight[0] = ct("0.08566224618958517252014807108636644676341125074202199119931771989947288027117007732396385271319433444");
      _points[1] = ct("0.9662347571015760139061507772469973045673828688561449124362748082633067504221000981381443699609629925");
      _weight[1] = ct("0.08566224618958517252014807108636644676341125074202199119931771989947288027117007732396385271319433508");
      _points[2] = ct("0.8306046932331322568306997975099526735032242821975850354072633529260917483035715504721432018732307282");
      _weight[2] = ct("0.1803807865240693037849167569188580558307609463733727411448696201185700189186308591604811009944096738");
      _points[3] = ct("0.1693953067668677431693002024900473264967757178024149645927366470739082516964284495278567981267692718");
      _weight[3] = ct("0.180380786524069303784916756918858055830760946373372741144869620118570018918630859160481100994409674");
      _points[4] = ct("0.6193095930415984543152508608403559677093053150700106750906975822871374671378199211246122136286745658");
      _weight[4] = ct("0.2339569672863455236949351719947754974058278028846052676558126599819571008101990635155550462923959912");
      _points[5] = ct("0.3806904069584015456847491391596440322906946849299893249093024177128625328621800788753877863713254342");
      _weight[5] = ct("0.2339569672863455236949351719947754974058278028846052676558126599819571008101990635155550462923959912");
      break;

    default :
      DUNE_THROW(QuadratureOrderOutOfRange, "Quadrature rule " << p << " not supported!");
    }
  }

} // namespace
