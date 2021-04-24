use crate::f64x2;

impl f64x2 {
    // |x| < pi / 4
    #[inline]
    pub(crate) fn sin_remez(self) -> Self {
        const C1: f64x2 = f64x2 { hi: 1.0, lo: -3.710932568160161e-33 };
        const C3: f64x2 = f64x2 { hi: -0.16666666666666666, lo: -9.251858538542761e-18 };
        const C5: f64x2 = f64x2 { hi: 0.008333333333333333, lo: 1.1564823172823584e-19 };
        const C7: f64x2 = f64x2 { hi: -0.0001984126984126984, lo: -1.7209555496532568e-22 };
        const C9: f64x2 = f64x2 { hi: 2.7557319223985893e-6, lo: -1.858394521208069e-22 };
        const C11: f64x2 = f64x2 { hi: -2.505210838544172e-8, lo: 1.4491633801389212e-24 };
        const C13: f64x2 = f64x2 { hi: 1.6059043836821613e-10, lo: 1.1931711466932894e-26 };
        const C15: f64x2 = f64x2 { hi: -7.647163731819808e-13, lo: 3.3092927077545036e-29 };
        const C17: f64x2 = f64x2 { hi: 2.811457254344742e-15, lo: 1.7646045165973752e-31 };
        const C19: f64x2 = f64x2 { hi: -8.22063524611474e-18, lo: -2.1399979358343384e-34 };
        const C21: f64x2 = f64x2 { hi: 1.957294082718338e-20, lo: -1.4638265350100345e-36 };
        const C23: f64x2 = f64x2 { hi: -3.868162594633509e-23, lo: 2.75299830507004e-39 };
        const C25: f64x2 = f64x2 { hi: 6.445350774747155e-26, lo: 3.53578137487626e-42 };
        const C27: f64x2 = f64x2 { hi: -8.983511414822678e-29, lo: -2.7266822373846815e-45 };

        let x = self;
        let x2 = x * x;
        let x4 = x2 * x2;
        let x6 = x2 * x4;
        let x8 = x4 * x4;
        let x16 = x8 * x8;
        let x24 = x8 * x16;

        let r1 = C1 + x2 * C3 + x4 * C5 + x6 * C7;
        let r2 = x8 * (C9 + x2 * C11 + x4 * C13 + x6 * C15);
        let r3 = x16 * (C17 + x2 * C19 + x4 * C21 + x6 * C23);
        let r4 = x24 * (C25 + x2 * C27);

        x * (r1 + r2 + r3 + r4)
    }

    // |x| < pi / 4
    #[inline]
    pub(crate) fn cos_remez(self) -> Self {
        const C0: f64x2 = f64x2 { hi: 1.0, lo: -7.415920157425322e-33 };
        const C2: f64x2 = f64x2 { hi: -0.5, lo: 1.178418360619204e-30 };
        const C4: f64x2 = f64x2 { hi: 0.041666666666666664, lo: 2.312964634604693e-18 };
        const C6: f64x2 = f64x2 { hi: -0.001388888888888889, lo: 5.300543986595874e-20 };
        const C8: f64x2 = f64x2 { hi: 2.48015873015873e-5, lo: 2.151020312414919e-23 };
        const C10: f64x2 = f64x2 { hi: -2.755731922398589e-7, lo: -2.3762056217062672e-23 };
        const C12: f64x2 = f64x2 { hi: 2.08767569878681e-9, lo: -1.3262132406487786e-25 };
        const C14: f64x2 = f64x2 { hi: -1.1470745597729708e-11, lo: 5.843851055946107e-28 };
        const C16: f64x2 = f64x2 { hi: 4.779477332385702e-14, lo: -3.0710780807342914e-30 };
        const C18: f64x2 = f64x2 { hi: -1.561920696740829e-16, lo: 8.75341824575357e-33 };
        const C20: f64x2 = f64x2 { hi: 4.1103175654770674e-19, lo: -6.63263393917438e-36 };
        const C22: f64x2 = f64x2 { hi: -8.89677188795462e-22, lo: -1.2759592625590192e-38 };
        const C24: f64x2 = f64x2 { hi: 1.6113071539190217e-24, lo: -2.9459654925058426e-41 };
        const C26: f64x2 = f64x2 { hi: -2.4235676574866644e-27, lo: 2.2034264640019724e-44 };

        let x = self;
        let x2 = x * x;
        let x4 = x2 * x2;
        let x6 = x2 * x4;
        let x8 = x4 * x4;
        let x16 = x8 * x8;
        let x24 = x8 * x16;

        let r1 = C0 + x2 * C2 + x4 * C4 + x6 * C6;
        let r2 = x8 * (C8 + x2 * C10 + x4 * C12 + x6 * C14);
        let r3 = x16 * (C16 + x2 * C18 + x4 * C20 + x6 * C22);
        let r4 = x24 * (C24 + x2 * C26);

        r1 + r2 + r3 + r4
    }

    // 1 / sqrt(2) < self < sqrt(2)
    #[inline]
    pub(crate) fn ln_remez(self) -> Self {
        assert!(0.707 <= self.hi && self.hi <= 1.415);
        const C1: f64x2 = f64x2 { hi: 2.0, lo: 2.531693403050348e-32 };
        const C3: f64x2 = f64x2 { hi: 0.6666666666666666, lo: 3.700743415403453e-17 };
        const C5: f64x2 = f64x2 { hi: 0.4, lo: -2.220446027080773e-17 };
        const C7: f64x2 = f64x2 { hi: 0.2857142857142857, lo: 1.586016141254861e-17 };
        const C9: f64x2 = f64x2 { hi: 0.2222222222222222, lo: 1.2407739663861083e-17 };
        const C11: f64x2 = f64x2 { hi: 0.1818181818181818, lo: 3.2065997154550064e-18 };
        const C13: f64x2 = f64x2 { hi: 0.1538461538461574, lo: -3.1965060154741107e-18 };
        const C15: f64x2 = f64x2 { hi: 0.13333333333287886, lo: 7.44820907006051e-18 };
        const C17: f64x2 = f64x2 { hi: 0.11764705886515148, lo: -1.17519451979169e-18 };
        const C19: f64x2 = f64x2 { hi: 0.1052631551291071, lo: 4.4806067994093946e-18 };
        const C21: f64x2 = f64x2 { hi: 0.095238228662164, lo: -3.577683177336467e-18 };
        const C23: f64x2 = f64x2 { hi: 0.08695190137972356, lo: -5.409228774321587e-18 };
        const C25: f64x2 = f64x2 { hi: 0.08011166789102804, lo: -2.092932986143882e-18 };
        const C27: f64x2 = f64x2 { hi: 0.07229374987314603, lo: -4.021248728080285e-18 };
        const C29: f64x2 = f64x2 { hi: 0.0855816562700506, lo: 6.675796100148133e-19 };

        let s = self;
        let x = (s - 1.0) / (s + 1.0);
        let x2 = x * x;
        let x4 = x2 * x2;
        let x6 = x2 * x4;
        let x8 = x4 * x4;
        let x10 = x4 * x6;
        let x20 = x10 * x10;

        let r1 = C1 + x2 * C3 + x4 * C5 + x6 * C7 + x8 * C9;
        let r2 = x10 * (C11 + x2 * C13 + x4 * C15 + x6 * C17 + x8 * C19);
        let r3 = x20 * (C21 + x2 * C23 + x4 * C25 + x6 * C27 + x8 * C29);

        x * (r1 + r2 + r3)
    }

    // |x| <= 0.34657359028 = ln(2) / 2
    #[inline]
    pub(crate) fn exp_remez(self) -> Self {
        const C2: f64x2 = f64x2 { hi: 0.16666666666666666, lo: 9.251858538542447e-18 };
        const C4: f64x2 = f64x2 { hi: -0.002777777777777778, lo: 1.0601087929995308e-19 };
        const C6: f64x2 = f64x2 { hi: 6.613756613756614e-5, lo: -4.460173646997389e-21 };
        const C8: f64x2 = f64x2 { hi: -1.6534391534391535e-6, lo: 7.121962972677988e-23 };
        const C10: f64x2 = f64x2 { hi: 4.1753513975736114e-8, lo: 1.158547249215353e-24 };
        const C12: f64x2 = f64x2 { hi: -1.0568380277354652e-9, lo: 5.58404280005523e-26 };
        const C14: f64x2 = f64x2 { hi: 2.6765073029312422e-11, lo: 9.12829168899536e-28 };
        const C16: f64x2 = f64x2 { hi: -6.779357355296518e-13, lo: -1.936460166148485e-29 };
        const C18: f64x2 = f64x2 { hi: 1.71700970829093e-14, lo: -1.033795194722353e-30 };
        const C20: f64x2 = f64x2 { hi: -4.278007864569552e-16, lo: 1.5907971763038017e-32 };

        let x = self;
        let x2 = x * x;
        let x4 = x2 * x2;
        let x6 = x2 * x4;
        let x8 = x4 * x4;
        let x14 = x8 * x6;

        let r1 = x2 * C2 + x4 * C4 + x6 * C6;
        let r2 = x8 * (C8 + x2 * C10 + x4 * C12);
        let r3 = x14 * (C14 + x2 * C16 + x4 * C18 + x6 * C20);
        // r = 2.0 + r1 + r2 + r3;
        // r = x (exp(x) + 1) / (exp(x) - 1) => exp(x) = 1 + 2r / (r - x)
        let c = x - (r1 + r2 + r3);
        1.0 + x + x * c / (2.0 - c)
    }

    pub(crate) fn atan_remez(self) -> Self {
        const C1: f64x2 = f64x2 { hi: 1.0, lo: -6.28445248852659e-31 };
        const C3: f64x2 = f64x2 { hi: -0.3333333333333333, lo: -1.850371707453854e-17 };
        const C5: f64x2 = f64x2 { hi: 0.2, lo: -1.1102233318388352e-17 };
        const C7: f64x2 = f64x2 { hi: -0.14285714285714285, lo: -7.928429626678543e-18 };
        const C9: f64x2 = f64x2 { hi: 0.1111111111111111, lo: 5.610977749577722e-18 };
        const C11: f64x2 = f64x2 { hi: -0.0909090909090908, lo: 4.45334301730286e-18 };
        const C13: f64x2 = f64x2 { hi: 0.0769230769230615, lo: 2.370246042413161e-18 };
        const C15: f64x2 = f64x2 { hi: -0.06666666666519032, lo: -3.3900293469993338e-18 };
        const C17: f64x2 = f64x2 { hi: 0.058823529310302304, lo: 2.4275518212554156e-18 };
        const C19: f64x2 = f64x2 { hi: -0.05263157387756347, lo: 1.5237747127147057e-18 };
        const C21: f64x2 = f64x2 { hi: 0.047618863095637315, lo: -2.8390751379840543e-19 };
        const C23: f64x2 = f64x2 { hi: -0.0434734148857163, lo: -7.686856741213805e-19 };
        const C25: f64x2 = f64x2 { hi: 0.0399103569181399, lo: 2.326566503469594e-18 };
        const C27: f64x2 = f64x2 { hi: -0.035922805494120864, lo: -1.4742055459849575e-18 };
        const C29: f64x2 = f64x2 { hi: 0.025993659875133343, lo: 4.697290175572255e-19 };

        let x = self;
        let x2 = x * x;
        let x4 = x2 * x2;
        let x6 = x2 * x4;
        let x8 = x4 * x4;
        let x10 = x4 * x6;
        let x20 = x10 * x10;

        let r1 = C1 + x2 * C3 + x4 * C5 + x6 * C7 + x8 * C9;
        let r2 = x10 * (C11 + x2 * C13 + x4 * C15 + x6 * C17 + x8 * C19);
        let r3 = x20 * (C21 + x2 * C23 + x4 * C25 + x6 * C27 + x8 * C29);

        x * (r1 + r2 + r3)
    }

    pub(crate) fn erfc_1_5_remez(self) -> Self {
        let x = self / 2.0 - Self { hi: 1.5, lo: 0.0 };
        let mut y = f64x2 { hi: -4.173011522948619e-19, lo: 2.151773608244596e-35 };
        y = y * x + f64x2 { hi: 1.3171951275905166e-18, lo: -3.262784621764736e-35 };
        y = y * x + f64x2 { hi: 4.70920134237599e-19, lo: 4.695025163371922e-35 };
        y = y * x + f64x2 { hi: -1.0539748979166572e-18, lo: 6.136789583505206e-35 };
        y = y * x + f64x2 { hi: -2.1628646933979583e-17, lo: 5.649456981647017e-34 };
        y = y * x + f64x2 { hi: 6.580080221557177e-17, lo: 6.332951861477313e-35 };
        y = y * x + f64x2 { hi: -1.2373823575923855e-16, lo: -1.164767818617041e-32 };
        y = y * x + f64x2 { hi: 3.912094378916367e-16, lo: -2.1697293677976843e-32 };
        y = y * x + f64x2 { hi: -1.3890187807213416e-15, lo: -1.2953343107844087e-32 };
        y = y * x + f64x2 { hi: 4.2134424078740605e-15, lo: -2.2336367937095188e-31 };
        y = y * x + f64x2 { hi: -1.2403584510539087e-14, lo: 6.425382951448193e-31 };
        y = y * x + f64x2 { hi: 3.712538218532184e-14, lo: -1.2162756920905071e-30 };
        y = y * x + f64x2 { hi: -1.104005200148566e-13, lo: -3.755768762339276e-30 };
        y = y * x + f64x2 { hi: 3.2392861490570504e-13, lo: 1.0714495435152694e-29 };
        y = y * x + f64x2 { hi: -9.40836582867467e-13, lo: -9.299673081208348e-29 };
        y = y * x + f64x2 { hi: 2.706443415804013e-12, lo: 6.912252009912485e-29 };
        y = y * x + f64x2 { hi: -7.705780443655945e-12, lo: 8.034132623545911e-28 };
        y = y * x + f64x2 { hi: 2.1708234830309646e-11, lo: 6.271618765709256e-28 };
        y = y * x + f64x2 { hi: -6.049556267480145e-11, lo: 1.1404323556150104e-27 };
        y = y * x + f64x2 { hi: 1.6672192180661353e-10, lo: -2.8401099739548477e-27 };
        y = y * x + f64x2 { hi: -4.5425553449845657e-10, lo: 1.590355659050682e-26 };
        y = y * x + f64x2 { hi: 1.223229666335777e-9, lo: 3.6179759042195166e-26 };
        y = y * x + f64x2 { hi: -3.2543929909606733e-9, lo: 3.062747000758059e-26 };
        y = y * x + f64x2 { hi: 8.551278439313862e-9, lo: 5.278287334501506e-25 };
        y = y * x + f64x2 { hi: -2.2183297526170795e-8, lo: -5.813643036256803e-25 };
        y = y * x + f64x2 { hi: 5.679096201163422e-8, lo: -1.5865659588563612e-25 };
        y = y * x + f64x2 { hi: -1.434175990187602e-7, lo: -6.462900281325292e-24 };
        y = y * x + f64x2 { hi: 3.571038035537446e-7, lo: 2.7411710850602287e-24 };
        y = y * x + f64x2 { hi: -8.762725030012075e-7, lo: -2.4107502745616608e-23 };
        y = y * x + f64x2 { hi: 2.117892312498399e-6, lo: -3.5360786953769844e-23 };
        y = y * x + f64x2 { hi: -5.038917537624995e-6, lo: 3.889677811095427e-22 };
        y = y * x + f64x2 { hi: 1.1794160931434194e-5, lo: -6.272867101557647e-22 };
        y = y * x + f64x2 { hi: -2.7139211780198174e-5, lo: -9.127206866079036e-22 };
        y = y * x + f64x2 { hi: 6.134859930030712e-5, lo: -5.301416563020298e-21 };
        y = y * x + f64x2 { hi: -0.0001361241180932827, lo: 1.0964761574067474e-20 };
        y = y * x + f64x2 { hi: 0.0002962090760903847, lo: -4.8995683948138694e-21 };
        y = y * x + f64x2 { hi: -0.0006314842765138408, lo: -1.8359426383544857e-20 };
        y = y * x + f64x2 { hi: 0.001317487759883742, lo: -5.645323949019005e-21 };
        y = y * x + f64x2 { hi: -0.002686651450903684, lo: 1.5044954873951244e-19 };
        y = y * x + f64x2 { hi: 0.005347464936239268, lo: -1.4480786246822472e-20 };
        y = y * x + f64x2 { hi: -0.010372017423899626, lo: 5.870452532873519e-19 };
        y = y * x + f64x2 { hi: 0.01956862483802889, lo: 4.0961416898882186e-19 };
        y = y * x + f64x2 { hi: -0.0358354481469806, lo: -3.067753354032795e-20 };
        y = y * x + f64x2 { hi: 0.06353748463948534, lo: 2.508233845053262e-19 };
        y = y * x + f64x2 { hi: -0.10874452001434574, lo: 1.3469843241290582e-18 };
        y = y * x + f64x2 { hi: 0.17900115118138996, lo: -5.4272175920200235e-18 };

        y
    }

    pub(crate) fn erfc_5_10_remez(self) -> Self {
        let x = self / 2.5 - Self { hi: 3.0, lo: 0.0 };
        let mut y = f64x2 { hi: 4.2991680017481634e-21, lo: -2.4989428780159913e-37 };
        y = y * x + f64x2 { hi: -1.5506536907573886e-20, lo: 1.4125641967754499e-36 };
        y = y * x + f64x2 { hi: 1.6890931272159943e-20, lo: -1.622347293480248e-37 };
        y = y * x + f64x2 { hi: -6.631216655043095e-20, lo: -4.058605045098766e-36 };
        y = y * x + f64x2 { hi: 4.170672248411603e-19, lo: 2.2868525293223118e-35 };
        y = y * x + f64x2 { hi: -1.497703484646083e-18, lo: 3.6559566676058296e-35 };
        y = y * x + f64x2 { hi: 4.950711200646159e-18, lo: 1.2906837301278596e-34 };
        y = y * x + f64x2 { hi: -1.7832129610866456e-17, lo: -2.356925974901086e-34 };
        y = y * x + f64x2 { hi: 6.451416952293557e-17, lo: 4.6394862044827265e-33 };
        y = y * x + f64x2 { hi: -2.2891991097758103e-16, lo: -2.2341667261894608e-32 };
        y = y * x + f64x2 { hi: 8.06457752377095e-16, lo: -4.445183120516124e-32 };
        y = y * x + f64x2 { hi: -2.826875054998268e-15, lo: -1.0389542680156556e-32 };
        y = y * x + f64x2 { hi: 9.845070747071263e-15, lo: -5.7729334938979695e-31 };
        y = y * x + f64x2 { hi: -3.405349886558733e-14, lo: -7.078413689058958e-31 };
        y = y * x + f64x2 { hi: 1.1697578418959786e-13, lo: -1.1068336248264454e-29 };
        y = y * x + f64x2 { hi: -3.989634711372613e-13, lo: 1.0343014251368376e-29 };
        y = y * x + f64x2 { hi: 1.3507295540717677e-12, lo: -2.1582211276480603e-30 };
        y = y * x + f64x2 { hi: -4.538257080442283e-12, lo: 2.884589994022212e-28 };
        y = y * x + f64x2 { hi: 1.512761587909949e-11, lo: -9.627112377046782e-28 };
        y = y * x + f64x2 { hi: -5.001185594604078e-11, lo: -3.958250640666782e-28 };
        y = y * x + f64x2 { hi: 1.6392272594079077e-10, lo: -1.0059690486863614e-26 };
        y = y * x + f64x2 { hi: -5.324638476070872e-10, lo: 6.565587905053832e-27 };
        y = y * x + f64x2 { hi: 1.7132484236062422e-9, lo: -9.837330468136443e-26 };
        y = y * x + f64x2 { hi: -5.4574678085162065e-9, lo: 2.469450519664579e-25 };
        y = y * x + f64x2 { hi: 1.719974809828779e-8, lo: 1.6499217118497196e-24 };
        y = y * x + f64x2 { hi: -5.358935050894508e-8, lo: -2.5588999529610755e-24 };
        y = y * x + f64x2 { hi: 1.649130365068191e-7, lo: -2.856162784012182e-24 };
        y = y * x + f64x2 { hi: -5.006694788416241e-7, lo: -4.568588521234318e-23 };
        y = y * x + f64x2 { hi: 1.4973926442280666e-6, lo: -8.990597381164192e-23 };
        y = y * x + f64x2 { hi: -4.403442164181749e-6, lo: -2.253398325570919e-22 };
        y = y * x + f64x2 { hi: 1.270084373012396e-5, lo: -2.9528515202596475e-22 };
        y = y * x + f64x2 { hi: -3.5805642403517175e-5, lo: 4.71428498204333e-22 };
        y = y * x + f64x2 { hi: 9.816805688521794e-5, lo: 5.96449459414328e-21 };
        y = y * x + f64x2 { hi: -0.0002597423991744517, lo: -2.42094445068534e-20 };
        y = y * x + f64x2 { hi: 0.0006547955362219539, lo: 3.83272594542379e-20 };
        y = y * x + f64x2 { hi: -0.0015356281323012648, lo: 5.3236490788364686e-20 };
        y = y * x + f64x2 { hi: 0.0031760235503980423, lo: -1.3503508956845972e-19 };
        y = y * x + f64x2 { hi: 0.5593026979715752, lo: -4.643965067263548e-17 };

        y
    }

    pub(crate) fn erfc_1_remez(self) -> Self {
        let x = self;
        let z = x * x;
        let mut y = f64x2 { hi: 1.4654075661333708e-19, lo: 9.860968581817586e-36 };
        y = y * z + f64x2 { hi: -4.3655412067160494e-18, lo: -1.8013039670415358e-34 };
        y = y * z + f64x2 { hi: 8.952761973695784e-17, lo: -4.8816915832622574e-33 };
        y = y * z + f64x2 { hi: -1.632084397548711e-15, lo: 8.964152287214925e-32 };
        y = y * z + f64x2 { hi: 2.7832025232513295e-14, lo: -7.114627633321976e-31 };
        y = y * z + f64x2 { hi: -4.463189973444574e-13, lo: 1.6294364728832262e-29 };
        y = y * z + f64x2 { hi: 6.711363958469121e-12, lo: -2.3273116044213122e-28 };
        y = y * z + f64x2 { hi: -9.422758873539076e-11, lo: 2.4224005951127012e-27 };
        y = y * z + f64x2 { hi: 1.2290555291825243e-9, lo: 6.683714882011236e-26 };
        y = y * z + f64x2 { hi: -1.4807192815477597e-8, lo: -1.798884593580444e-25 };
        y = y * z + f64x2 { hi: 1.63658446912222e-7, lo: 8.091381745719929e-24 };
        y = y * z + f64x2 { hi: -1.6462114365888935e-6, lo: -9.425224149920854e-23 };
        y = y * z + f64x2 { hi: 1.4925650358406245e-5, lo: -9.082422345266634e-23 };
        y = y * z + f64x2 { hi: -0.00012055332981789664, lo: -5.677432254336662e-21 };
        y = y * z + f64x2 { hi: 0.0008548327023450853, lo: -5.022895506135341e-20 };
        y = y * z + f64x2 { hi: -0.005223977625442188, lo: 8.963053564118845e-20 };
        y = y * z + f64x2 { hi: 0.026866170645131252, lo: -4.60929048899779e-19 };
        y = y * z + f64x2 { hi: -0.11283791670955126, lo: 4.01756916792255e-18 };
        y = y * z + f64x2 { hi: 0.37612638903183754, lo: -1.3391897206103609e-17 };
        y = y * z + f64x2 { hi: -1.1283791670955126, lo: -1.533545961316562e-17 };
        y = y * x;

        y
    }
}
