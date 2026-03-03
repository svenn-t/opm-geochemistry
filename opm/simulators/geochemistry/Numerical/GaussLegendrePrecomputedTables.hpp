/*
 * MIT License
 *
 * Copyright (C) 2025 Aksel Hiorth
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/
#ifndef GAUSS_LEGENDRE_PRECOMPUTED_TABLES_IS_DEF_H
#define GAUSS_LEGENDRE_PRECOMPUTED_TABLES_IS_DEF_H

#include <opm/simulators/geochemistry/Common/CustomExceptions.hpp>

#include <array>

/*
* Hard-coded values are taken from a C++ project made by Pavel Holoborodko:
*       http://www.holoborodko.com/pavel/numerical-methods/numerical-integration/.
*/
template <class Real, unsigned N>
class GaussLegendreTables{
public:
    // TO DO(maybe): Implement general case by computing the needed numbers on-the-fly...?
    static const std::vector<Real>& abscissas(){
        throw NotImplementedException("Only hard-coded versions of the GaussLegendre integration method can be used!");

    }
    static const std::vector<Real>& weights()
    {
        throw NotImplementedException("Only hard-coded versions of the GaussLegendre integration method can be used!");
    }
};

// ====================================== EVEN ORDER CASES ======================================
template <class T>
class GaussLegendreTables<T, 2>{
public:
    static std::array<T, 1> const& abscissas(){
        static std::array<T, 1> data = { 0.5773502691896257645091488 };
        return data;
    }
    static std::array<T, 1> const& weights(){
        static std::array<T, 1> data = { 1.0000000000000000000000000 };
        return data;
    }
};

template <class T>
class GaussLegendreTables<T, 4>{
public:
    static std::array<T, 2> const& abscissas(){
        static std::array<T, 2> data = { 0.3399810435848562648026658, 0.8611363115940525752239465 };
        return data;
    }
    static std::array<T, 2> const& weights(){
        static std::array<T, 2> data = { 0.6521451548625461426269361, 0.3478548451374538573730639 };
        return data;
    }
};

template <class T>
class GaussLegendreTables<T, 6>{
public:
    static std::array<T, 3> const& abscissas(){
        static std::array<T, 3> data = { 0.2386191860831969086305017, 0.6612093864662645136613996, 0.9324695142031520278123016 };
        return data;
    }
    static std::array<T, 3> const& weights(){
        static std::array<T, 3> data = { 0.4679139345726910473898703, 0.3607615730481386075698335, 0.1713244923791703450402961 };
        return data;
    }
};

template <class T>
class GaussLegendreTables<T, 8>{
public:
    static std::array<T, 4> const& abscissas(){
        static std::array<T, 4> data = { 0.1834346424956498049394761, 0.5255324099163289858177390, 0.7966664774136267395915539, 0.9602898564975362316835609 };
        return data;
    }
    static std::array<T, 4> const& weights(){
        static std::array<T, 4> data = { 0.3626837833783619829651504, 0.3137066458778872873379622, 0.2223810344533744705443560, 0.1012285362903762591525314 };
        return data;
    }
};

template <class T>
class GaussLegendreTables<T, 10>{
public:
    static std::array<T, 5> const& abscissas(){
        static std::array<T, 5> data = { 0.1488743389816312108848260, 0.4333953941292471907992659, 0.6794095682990244062343274, 0.8650633666889845107320967, 0.9739065285171717200779640 };
        return data;
    }
    static std::array<T, 5> const& weights(){
        static std::array<T, 5> data = { 0.2955242247147528701738930, 0.2692667193099963550912269, 0.2190863625159820439955349, 0.1494513491505805931457763, 0.0666713443086881375935688 };
        return data;
    }
};

template <class T>
class GaussLegendreTables<T, 12>{
public:
    static std::array<T, 6> const& abscissas(){
        static std::array<T, 6> data =
        {
            0.1252334085114689154724414,
            0.3678314989981801937526915,
            0.5873179542866174472967024,
            0.7699026741943046870368938,
            0.9041172563704748566784659,
            0.9815606342467192506905491
        };
        return data;
    }
    static std::array<T, 6> const& weights(){
        static std::array<T, 6> data =
        {
            0.2491470458134027850005624,
            0.2334925365383548087608499,
            0.2031674267230659217490645,
            0.1600783285433462263346525,
            0.1069393259953184309602547,
            0.0471753363865118271946160
        };
        return data;
    }
};

template <class T>
class GaussLegendreTables<T, 14>{
public:
    static std::array<T, 7> const& abscissas(){
        static std::array<T, 7> data =
        {
            0.1080549487073436620662447,
            0.3191123689278897604356718,
            0.5152486363581540919652907,
            0.6872929048116854701480198,
            0.8272013150697649931897947,
            0.9284348836635735173363911,
            0.9862838086968123388415973
        };
        return data;
    }
    static std::array<T, 7> const& weights(){
        static std::array<T, 7> data =
        {
            0.2152638534631577901958764,
            0.2051984637212956039659241,
            0.1855383974779378137417166,
            0.1572031671581935345696019,
            0.1215185706879031846894148,
            0.0801580871597602098056333,
            0.0351194603317518630318329
        };
        return data;
    }
};

template <class T>
class GaussLegendreTables<T, 16>{
public:
    static std::array<T, 8> const& abscissas(){
        static std::array<T, 8> data =
        {
            0.0950125098376374401853193,
            0.2816035507792589132304605,
            0.4580167776572273863424194,
            0.6178762444026437484466718,
            0.7554044083550030338951012,
            0.8656312023878317438804679,
            0.9445750230732325760779884,
            0.9894009349916499325961542
        };
        return data;
    }
    static std::array<T, 8> const& weights(){
        static std::array<T, 8> data =
        {
            0.1894506104550684962853967,
            0.1826034150449235888667637,
            0.1691565193950025381893121,
            0.1495959888165767320815017,
            0.1246289712555338720524763,
            0.0951585116824927848099251,
            0.0622535239386478928628438,
            0.0271524594117540948517806
        };
        return data;
    }
};

template <class T>
class GaussLegendreTables<T, 18>{
public:
    static std::array<T, 9> const& abscissas(){
        static std::array<T, 9> data =
        {
            0.0847750130417353012422619,
            0.2518862256915055095889729,
            0.4117511614628426460359318,
            0.5597708310739475346078715,
            0.6916870430603532078748911,
            0.8037049589725231156824175,
            0.8926024664975557392060606,
            0.9558239495713977551811959,
            0.9915651684209309467300160
        };
        return data;
    }
    static std::array<T, 9> const& weights(){
        static std::array<T, 9> data =
        {
            0.1691423829631435918406565,
            0.1642764837458327229860538,
            0.1546846751262652449254180,
            0.1406429146706506512047313,
            0.1225552067114784601845191,
            0.1009420441062871655628140,
            0.0764257302548890565291297,
            0.0497145488949697964533349,
            0.0216160135264833103133427
        };
        return data;
    }
};

template <class T>
class GaussLegendreTables<T, 20>{
public:
    static std::array<T, 10> const& abscissas(){
        static std::array<T, 10> data =
        {
            0.0765265211334973337546404,
            0.2277858511416450780804962,
            0.3737060887154195606725482,
            0.5108670019508270980043641,
            0.6360536807265150254528367,
            0.7463319064601507926143051,
            0.8391169718222188233945291,
            0.9122344282513259058677524,
            0.9639719272779137912676661,
            0.9931285991850949247861224
        };
        return data;
    }
    static std::array<T, 10> const& weights(){
        static std::array<T, 10> data =
        {
            0.1527533871307258506980843,
            0.1491729864726037467878287,
            0.1420961093183820513292983,
            0.1316886384491766268984945,
            0.1181945319615184173123774,
            0.1019301198172404350367501,
            0.0832767415767047487247581,
            0.0626720483341090635695065,
            0.0406014298003869413310400,
            0.0176140071391521183118620
        };
        return data;
    }
};

// ======================================= ODD ORDER CASES ======================================

template <class T>
class GaussLegendreTables<T, 3>{
public:
    static std::array<T, 2> const& abscissas(){
        static std::array<T, 2> data = { 0.0000000000000000000000000, 0.7745966692414833770358531 };
        return data;
    }
    static std::array<T, 2> const& weights(){
        static std::array<T, 2> data = { 0.8888888888888888888888889, 0.5555555555555555555555556 };
        return data;
    }
};

template <class T>
class GaussLegendreTables<T, 5>{
public:
    static std::array<T, 3> const& abscissas(){
        static std::array<T, 3> data = { 0.0000000000000000000000000, 0.5384693101056830910363144, 0.9061798459386639927976269 };
        return data;
    }
    static std::array<T, 3> const& weights(){
        static std::array<T, 3> data = { 0.5688888888888888888888889, 0.4786286704993664680412915, 0.2369268850561890875142640 };
        return data;
    }
};


template <class T>
class GaussLegendreTables<T, 7>{
public:
    static std::array<T, 4> const& abscissas(){
        static std::array<T, 4> data = { 0.0000000000000000000000000, 0.4058451513773971669066064, 0.7415311855993944398638648, 0.9491079123427585245261897 };
        return data;
    }
    static std::array<T, 4> const& weights(){
        static std::array<T, 4> data = { 0.4179591836734693877551020, 0.3818300505051189449503698, 0.2797053914892766679014678, 0.1294849661688696932706114 };
        return data;
    }
};

template <class T>
class GaussLegendreTables<T, 9>{
public:
    static std::array<T, 5> const& abscissas(){
        static std::array<T, 5> data =
        {
            0.0000000000000000000000000,
            0.3242534234038089290385380,
            0.6133714327005903973087020,
            0.8360311073266357942994298,
            0.9681602395076260898355762
        };
        return data;
    }
    static std::array<T, 5> const& weights(){
        static std::array<T, 5> data =
        {
            0.3302393550012597631645251,
            0.3123470770400028400686304,
            0.2606106964029354623187429,
            0.1806481606948574040584720,
            0.0812743883615744119718922
        };
        return data;
    }
};

template <class T>
class GaussLegendreTables<T, 11>{
public:
    static std::array<T, 6> const& abscissas(){
        static std::array<T, 6> data =
        {
        0.0000000000000000000000000,
        0.2695431559523449723315320,
        0.5190961292068118159257257,
        0.7301520055740493240934163,
        0.8870625997680952990751578,
        0.9782286581460569928039380
        };
        return data;
    }
    static std::array<T, 6> const& weights(){
        static std::array<T, 6> data =
        {
            0.2729250867779006307144835,
            0.2628045445102466621806889,
            0.2331937645919904799185237,
            0.1862902109277342514260976,
            0.1255803694649046246346943,
            0.0556685671161736664827537
        };
        return data;
    }
};

template <class T>
class GaussLegendreTables<T, 13>{
public:
    static std::array<T, 7> const& abscissas(){
        static std::array<T, 7> data =
        {
            0.0000000000000000000000000,
            0.2304583159551347940655281,
            0.4484927510364468528779129,
            0.6423493394403402206439846,
            0.8015780907333099127942065,
            0.9175983992229779652065478,
            0.9841830547185881494728294
        };
        return data;
    }
    static std::array<T, 7> const& weights(){
        static std::array<T, 7> data =
        {
            0.2325515532308739101945895,
            0.2262831802628972384120902,
            0.2078160475368885023125232,
            0.1781459807619457382800467,
            0.1388735102197872384636018,
            0.0921214998377284479144218,
            0.0404840047653158795200216
        };
        return data;
    }
};

template <class T>
class GaussLegendreTables<T, 15>{
public:
    static std::array<T, 8> const& abscissas(){
        static std::array<T, 8> data =
        {
            0.0000000000000000000000000,
            0.2011940939974345223006283,
            0.3941513470775633698972074,
            0.5709721726085388475372267,
            0.7244177313601700474161861,
            0.8482065834104272162006483,
            0.9372733924007059043077589,
            0.9879925180204854284895657
        };
        return data;
    }
    static std::array<T, 8> const& weights(){
        static std::array<T, 8> data =
        {
            0.2025782419255612728806202,
            0.1984314853271115764561183,
            0.1861610000155622110268006,
            0.1662692058169939335532009,
            0.1395706779261543144478048,
            0.1071592204671719350118695,
            0.0703660474881081247092674,
            0.0307532419961172683546284
        };
        return data;
    }
};

template <class T>
class GaussLegendreTables<T, 17>{
public:
    static std::array<T, 9> const& abscissas(){
        static std::array<T, 9> data =
        {
            0.0000000000000000000000000,
            0.1784841814958478558506775,
            0.3512317634538763152971855,
            0.5126905370864769678862466,
            0.6576711592166907658503022,
            0.7815140038968014069252301,
            0.8802391537269859021229557,
            0.9506755217687677612227170,
            0.9905754753144173356754340
        };
        return data;
    }
    static std::array<T, 9> const& weights(){
        static std::array<T, 9> data =
        {
            0.1794464703562065254582656,
            0.1765627053669926463252710,
            0.1680041021564500445099707,
            0.1540457610768102880814316,
            0.1351363684685254732863200,
            0.1118838471934039710947884,
            0.0850361483171791808835354,
            0.0554595293739872011294402,
            0.0241483028685479319601100
        };
        return data;
    }
};

template <class T>
class GaussLegendreTables<T, 19>{
public:
    static std::array<T, 10> const& abscissas(){
        static std::array<T, 10> data =
        {
            0.0000000000000000000000000,
            0.1603586456402253758680961,
            0.3165640999636298319901173,
            0.4645707413759609457172671,
            0.6005453046616810234696382,
            0.7209661773352293786170959,
            0.8227146565371428249789225,
            0.9031559036148179016426609,
            0.9602081521348300308527788,
            0.9924068438435844031890177
        };
        return data;
    }
    static std::array<T, 10> const& weights(){
        static std::array<T, 10> data =
        {
            0.1610544498487836959791636,
            0.1589688433939543476499564,
            0.1527660420658596667788554,
            0.1426067021736066117757461,
            0.1287539625393362276755158,
            0.1115666455473339947160239,
            0.0914900216224499994644621,
            0.0690445427376412265807083,
            0.0448142267656996003328382,
            0.0194617882297264770363120
        };
        return data;
    }
};

#endif
