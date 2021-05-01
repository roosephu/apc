use num::Complex;

use F64x2::f64x2;

pub(crate) const RS_GABCKE_GAMMA: [[f64x2; 26]; 11] = [
    [
        f64x2 { hi: 0.6426672862397684, lo: -5.405020898276905e-17 },
        f64x2 { hi: 0.27197299999785507, lo: -2.324922478643637e-18 },
        f64x2 { hi: 0.010738605819340285, lo: -7.259826353252223e-19 },
        f64x2 { hi: -0.0013743815296336614, lo: -3.989397266596713e-20 },
        f64x2 { hi: -0.00012468221880320676, lo: -8.116850720037552e-21 },
        f64x2 { hi: -5.764599706783048e-7, lo: -1.3497354112417146e-23 },
        f64x2 { hi: 2.728067429580452e-7, lo: 2.531641491695852e-23 },
        f64x2 { hi: 8.07795305950047e-9, lo: -1.6348465461165598e-25 },
        f64x2 { hi: -2.0884608068869654e-10, lo: -9.186698446291329e-27 },
        f64x2 { hi: -1.3115561854739528e-11, lo: 7.26645367460136e-28 },
        f64x2 { hi: -1.4207987228087186e-14, lo: 3.5613660611352426e-31 },
        f64x2 { hi: 1.0271701357931162e-14, lo: -1.0581401316172562e-32 },
        f64x2 { hi: 1.3974598819518373e-16, lo: 1.193118980994164e-32 },
        f64x2 { hi: -4.4841187339522885e-18, lo: 1.5721523238671847e-34 },
        f64x2 { hi: -1.1830599573845289e-19, lo: -3.9856164969516115e-36 },
        f64x2 { hi: 9.389869560399935e-22, lo: 3.281591204345519e-38 },
        f64x2 { hi: 5.601822847320697e-23, lo: 6.8417995938899474e-40 },
        f64x2 { hi: 1.0023543875614807e-25, lo: 5.348697690815155e-42 },
        f64x2 { hi: -1.7592985581293886e-26, lo: -5.5994115608571294e-43 },
        f64x2 { hi: -1.4854553062733663e-28, lo: 2.7064838367273312e-45 },
        f64x2 { hi: 3.8087608010848314e-30, lo: -3.419287919722399e-46 },
        f64x2 { hi: 5.901183139529734e-32, lo: 3.748092600777806e-48 },
        f64x2 { hi: -5.364406588951707e-34, lo: -7.744542595962487e-51 },
        f64x2 { hi: -1.5079228966504238e-35, lo: -5.409955158384653e-52 },
        f64x2 { hi: 2.973747761357669e-38, lo: -1.5363277840772604e-54 },
        f64x2 { hi: 2.831060640801242e-39, lo: -1.3817188335567144e-55 },
    ],
    [
        f64x2 { hi: 0.003669156562182772, lo: 4.8517624939411e-20 },
        f64x2 { hi: -0.028734140966371547, lo: 1.5518517595356827e-18 },
        f64x2 { hi: -0.005607161520384222, lo: -5.625595236948406e-20 },
        f64x2 { hi: 2.0739220807279964e-5, lo: -4.0844886552168817e-22 },
        f64x2 { hi: 5.201208663127012e-5, lo: 1.5566761780892368e-21 },
        f64x2 { hi: 2.205823831031653e-6, lo: -1.692189267618361e-22 },
        f64x2 { hi: -1.0907385768109821e-7, lo: 2.6146175524798896e-24 },
        f64x2 { hi: -8.655485649453228e-9, lo: 1.0623903240733617e-25 },
        f64x2 { hi: 9.551112447669321e-12, lo: 4.6551394504185185e-28 },
        f64x2 { hi: 1.3188070728878103e-11, lo: -7.3484330408524975e-28 },
        f64x2 { hi: 2.115970918325531e-13, lo: 5.776893301730298e-30 },
        f64x2 { hi: -9.997138776383605e-15, lo: -4.857337723505849e-31 },
        f64x2 { hi: -3.078372420606269e-16, lo: 1.8649271059310332e-32 },
        f64x2 { hi: 3.4981526238874625e-18, lo: 3.5452582638741465e-34 },
        f64x2 { hi: 2.2574034284995735e-19, lo: -6.287600891713198e-36 },
        f64x2 { hi: 2.96580834967959e-22, lo: -5.4702528309363674e-39 },
        f64x2 { hi: -1.0357470543081694e-22, lo: -4.427876092094538e-39 },
        f64x2 { hi: -9.575082777353988e-25, lo: -6.684290595694097e-42 },
        f64x2 { hi: 3.1498466813118526e-26, lo: 5.864734462408796e-43 },
        f64x2 { hi: 5.377443845416828e-28, lo: -1.287342735441815e-45 },
        f64x2 { hi: -6.103388904262898e-30, lo: 3.0676106781300854e-46 },
        f64x2 { hi: -1.8459669844569243e-31, lo: 9.204996819578101e-49 },
        f64x2 { hi: 5.0179727598745494e-34, lo: -1.8095574950982408e-50 },
        f64x2 { hi: 4.519724243443113e-35, lo: -1.889030784568232e-51 },
        f64x2 { hi: 1.1716947888275927e-37, lo: 3.121823781193407e-54 },
        f64x2 { hi: -8.308858951569504e-39, lo: 1.2546811641908664e-55 },
    ],
    [
        f64x2 { hi: 0.0031461158539889122, lo: 1.2465025418958422e-20 },
        f64x2 { hi: -0.0023087838845307503, lo: 1.8692677027356787e-19 },
        f64x2 { hi: 5.769820766689844e-5, lo: -3.1495975520988944e-21 },
        f64x2 { hi: 0.000352388620236659, lo: 1.641690435167814e-20 },
        f64x2 { hi: 2.5246667458684434e-5, lo: 1.245610621231287e-22 },
        f64x2 { hi: -3.442821197193136e-6, lo: -7.366245642841296e-23 },
        f64x2 { hi: -3.535074556622459e-7, lo: 1.921491644506255e-23 },
        f64x2 { hi: 3.730830183792625e-9, lo: 2.1175449687662407e-25 },
        f64x2 { hi: 1.2776951864116635e-9, lo: -6.934192898379383e-27 },
        f64x2 { hi: 2.1874616204147057e-11, lo: 5.379313982478304e-28 },
        f64x2 { hi: -1.914141096461037e-12, lo: -5.427474616466113e-29 },
        f64x2 { hi: -6.562883102168523e-14, lo: 2.4190790934922127e-31 },
        f64x2 { hi: 1.2586009182411715e-15, lo: 4.825689145715241e-32 },
        f64x2 { hi: 8.140076623881463e-17, lo: -2.624278821206968e-33 },
        f64x2 { hi: -5.423874275488608e-20, lo: 3.5259410449687945e-36 },
        f64x2 { hi: -5.796980131086543e-20, lo: -4.0644353740097735e-36 },
        f64x2 { hi: -5.382916503746397e-22, lo: 6.802721571163337e-39 },
        f64x2 { hi: 2.6010080772383425e-23, lo: 1.3145995274166021e-39 },
        f64x2 { hi: 4.666966774911328e-25, lo: -4.509257845846034e-41 },
        f64x2 { hi: -7.288849536075177e-27, lo: -4.999257295701152e-43 },
        f64x2 { hi: -2.250096790723193e-28, lo: -1.3016407581648374e-44 },
        f64x2 { hi: 9.737854958612028e-31, lo: 6.810052849357968e-47 },
        f64x2 { hi: 7.414681256143464e-32, lo: 2.685750987629922e-48 },
        f64x2 { hi: 1.536964535218762e-34, lo: -3.8751975395748384e-51 },
        f64x2 { hi: -1.7848849203395843e-35, lo: 8.308815973775812e-52 },
        f64x2 { hi: -1.2518534861621486e-37, lo: -4.494787722826375e-54 },
    ],
    [
        f64x2 { hi: 0.00030191243459800856, lo: -2.605115738122299e-20 },
        f64x2 { hi: -0.0007462899936200945, lo: 5.8938721220285586e-21 },
        f64x2 { hi: 0.0002816038876567984, lo: -2.4794899748733012e-20 },
        f64x2 { hi: -2.3005646747348866e-5, lo: 5.415457062132067e-23 },
        f64x2 { hi: -1.3143346079994013e-5, lo: 3.8827689034409128e-22 },
        f64x2 { hi: 9.097570555313324e-8, lo: -1.5436008441642624e-24 },
        f64x2 { hi: 1.4295160201730648e-7, lo: 1.274537592500182e-23 },
        f64x2 { hi: 4.037920513056026e-9, lo: -4.698375276803934e-26 },
        f64x2 { hi: -4.877384974746115e-10, lo: 2.597849969291361e-26 },
        f64x2 { hi: -2.3372094790693514e-11, lo: 1.5974787174507963e-27 },
        f64x2 { hi: 6.188215896189135e-13, lo: 3.2398337539573857e-29 },
        f64x2 { hi: 5.1151190087141844e-14, lo: 1.5177375673377057e-31 },
        f64x2 { hi: -7.643137881406036e-17, lo: 1.306174308180276e-33 },
        f64x2 { hi: -5.889863661237705e-17, lo: 1.3801796593584514e-33 },
        f64x2 { hi: -6.3913282746240375e-19, lo: 3.688425627938971e-35 },
        f64x2 { hi: 4.008866571106544e-20, lo: -2.9773310250879572e-36 },
        f64x2 { hi: 8.337112840847493e-22, lo: -3.131613568675373e-38 },
        f64x2 { hi: -1.6325993418950848e-23, lo: 7.365791749696725e-40 },
        f64x2 { hi: -5.692307992769963e-25, lo: 1.9728906788374198e-41 },
        f64x2 { hi: 3.1769103118712025e-27, lo: -9.689478032574735e-44 },
        f64x2 { hi: 2.5560087721966874e-28, lo: 2.193229837379676e-45 },
        f64x2 { hi: 5.09370343454115e-31, lo: 5.0939035114023804e-48 },
        f64x2 { hi: -8.14437507813587e-32, lo: -4.3434353664568717e-48 },
        f64x2 { hi: -6.1846389085793075e-34, lo: 3.4255440621288768e-50 },
        f64x2 { hi: 1.8905839011132534e-35, lo: 1.0615317072505247e-51 },
        f64x2 { hi: 2.497523559493801e-37, lo: 1.9494099091736543e-53 },
    ],
    [
        f64x2 { hi: 0.0001676574524669686, lo: -1.09136640699526e-20 },
        f64x2 { hi: -0.00022728768943416726, lo: 3.14141358511196e-21 },
        f64x2 { hi: 6.477387188445696e-5, lo: -3.768685095371224e-21 },
        f64x2 { hi: -8.49220050012541e-6, lo: 4.1688523384080164e-22 },
        f64x2 { hi: -2.6161407245219076e-6, lo: -9.450416571045903e-23 },
        f64x2 { hi: 8.336764968733215e-7, lo: -4.2035464471276434e-23 },
        f64x2 { hi: 6.324704037544833e-8, lo: -3.7286766746352895e-24 },
        f64x2 { hi: -1.0059949403001072e-8, lo: 9.662938012848006e-26 },
        f64x2 { hi: -7.822677204130333e-10, lo: 3.774595892369452e-27 },
        f64x2 { hi: 3.16765828534986e-11, lo: 2.9694710589173263e-28 },
        f64x2 { hi: 3.5006944702052894e-12, lo: 5.738331548480427e-29 },
        f64x2 { hi: -1.4314814511443748e-14, lo: -1.1872704094602542e-30 },
        f64x2 { hi: -7.269402707921764e-15, lo: 2.714372772987968e-31 },
        f64x2 { hi: -8.780556594835957e-17, lo: 3.0214638264210787e-33 },
        f64x2 { hi: 8.15025447495458e-18, lo: -6.903631296516382e-34 },
        f64x2 { hi: 1.920839705822086e-19, lo: 6.326533872117065e-36 },
        f64x2 { hi: -5.175655213952982e-21, lo: 1.8494633231847035e-37 },
        f64x2 { hi: -1.976773672440578e-22, lo: -1.016067371237183e-38 },
        f64x2 { hi: 1.6059867343952904e-24, lo: -3.4527333411347024e-41 },
        f64x2 { hi: 1.2658632662439863e-25, lo: -7.612655457854657e-42 },
        f64x2 { hi: 1.6326189825047856e-28, lo: 9.18920417727519e-45 },
        f64x2 { hi: -5.537211021742346e-29, lo: 3.1115924096825166e-45 },
        f64x2 { hi: -4.310484397202635e-31, lo: 4.7775638185694345e-48 },
        f64x2 { hi: 1.7174510567185548e-32, lo: -2.494364832478607e-49 },
        f64x2 { hi: 2.384697343638105e-34, lo: 1.0812959357504413e-50 },
        f64x2 { hi: -3.7541877585640003e-36, lo: -2.0523485916007018e-53 },
    ],
    [
        f64x2 { hi: -0.0001009490716640513, lo: 1.9725490177462588e-21 },
        f64x2 { hi: 2.5321238631924578e-5, lo: -8.641975321309539e-24 },
        f64x2 { hi: 5.936131306732197e-6, lo: -6.698259415735041e-23 },
        f64x2 { hi: -5.569282352788995e-6, lo: -4.0712764197025003e-22 },
        f64x2 { hi: 1.3498287778014868e-6, lo: 8.218251367260146e-23 },
        f64x2 { hi: -1.842554298220938e-8, lo: -1.5380142042329987e-24 },
        f64x2 { hi: -3.700393942792748e-8, lo: -2.9135165026317216e-25 },
        f64x2 { hi: 7.814406763977287e-10, lo: -4.7930627201900524e-26 },
        f64x2 { hi: 3.7173748594546684e-10, lo: -2.1895945562246047e-26 },
        f64x2 { hi: 1.7631825761962343e-12, lo: -3.179665102012118e-29 },
        f64x2 { hi: -1.5421503978543737e-12, lo: -6.532338853339962e-29 },
        f64x2 { hi: -3.197827575699093e-14, lo: 1.0320424064514288e-30 },
        f64x2 { hi: 3.06157376568069e-15, lo: -1.5573289423494478e-31 },
        f64x2 { hi: 1.0134461604121639e-16, lo: 2.345659651129789e-33 },
        f64x2 { hi: -3.1318394339406034e-18, lo: 9.781960085837102e-36 },
        f64x2 { hi: -1.5700081019275403e-19, lo: -8.92863192704536e-36 },
        f64x2 { hi: 1.4404518422893965e-21, lo: -8.90346036503584e-38 },
        f64x2 { hi: 1.4599353364699049e-22, lo: 8.608203216122158e-39 },
        f64x2 { hi: 2.5960849210556536e-25, lo: 1.6745248521679822e-41 },
        f64x2 { hi: -8.929442752805088e-26, lo: -5.613450198578382e-42 },
        f64x2 { hi: -8.191704573158245e-28, lo: -2.925956986974212e-44 },
        f64x2 { hi: 3.7496284130881265e-29, lo: 2.7210175129674097e-45 },
        f64x2 { hi: 5.985960155147785e-31, lo: -3.2982098602754204e-47 },
        f64x2 { hi: -1.0810867382178206e-32, lo: -2.742692181719002e-49 },
        f64x2 { hi: -2.700212667743483e-34, lo: 3.726553200638024e-51 },
        f64x2 { hi: 1.950283166970115e-36, lo: 1.406083958881213e-52 },
    ],
    [
        f64x2 { hi: 1.2189742141068971e-5, lo: -5.2225751533143015e-22 },
        f64x2 { hi: -1.3829760140503787e-5, lo: 5.003786278007855e-22 },
        f64x2 { hi: 5.11096730499826e-6, lo: -3.3180918347775264e-22 },
        f64x2 { hi: -2.0458136450386076e-6, lo: -5.0893764664704725e-23 },
        f64x2 { hi: 4.938136644832012e-7, lo: -4.827516571898731e-23 },
        f64x2 { hi: -3.6187528349622816e-8, lo: 2.0460320389045594e-24 },
        f64x2 { hi: -1.287690509807986e-8, lo: -4.14760087085482e-25 },
        f64x2 { hi: 2.574412111144866e-9, lo: 1.3890852072700892e-25 },
        f64x2 { hi: 1.3641457070791684e-10, lo: -8.80575914978285e-27 },
        f64x2 { hi: -3.032439574084382e-11, lo: -1.4905970355761775e-27 },
        f64x2 { hi: -1.3216671239902537e-12, lo: 2.4181277347626273e-29 },
        f64x2 { hi: 1.3031652130009368e-13, lo: 3.8024965144209216e-30 },
        f64x2 { hi: 6.63588355320067e-15, lo: -1.0821927929883243e-31 },
        f64x2 { hi: -2.46003565479328e-16, lo: 6.82509636804775e-33 },
        f64x2 { hi: -1.6815279208168834e-17, lo: 2.798016668013323e-34 },
        f64x2 { hi: 1.8937932080359403e-19, lo: 5.303317249246921e-37 },
        f64x2 { hi: 2.430650612737236e-20, lo: -7.0032247615631054e-37 },
        f64x2 { hi: 4.608486141193199e-23, lo: 9.606969149142508e-41 },
        f64x2 { hi: -2.1956897626337115e-23, lo: -4.222398619067676e-41 },
        f64x2 { hi: -2.295880332596837e-25, lo: 5.4879181570433406e-42 },
        f64x2 { hi: 1.306599063382408e-26, lo: 9.343597045745804e-43 },
        f64x2 { hi: 2.3471644895116074e-28, lo: -1.472025868683126e-44 },
        f64x2 { hi: -5.181499667299091e-30, lo: 2.0890112735119453e-47 },
        f64x2 { hi: -1.4258922561649248e-31, lo: 4.092146894861322e-48 },
        f64x2 { hi: 1.2733209778170061e-33, lo: -5.035716915667804e-50 },
        f64x2 { hi: 6.029143930316113e-35, lo: 5.871277647519953e-52 },
    ],
    [
        f64x2 { hi: -1.827285786561044e-5, lo: -6.806781403601642e-23 },
        f64x2 { hi: 1.1008400136344443e-5, lo: 5.154690856384144e-22 },
        f64x2 { hi: -3.282532467061245e-6, lo: 1.7797999526499082e-22 },
        f64x2 { hi: 5.437662797676692e-7, lo: 2.1090476434528882e-23 },
        f64x2 { hi: 9.174553888200617e-9, lo: -7.018921323297588e-25 },
        f64x2 { hi: -2.9741427534891034e-8, lo: -1.5006719621650694e-24 },
        f64x2 { hi: 6.231294398552861e-9, lo: 4.532458130938043e-26 },
        f64x2 { hi: -1.2119656685887065e-10, lo: -2.5811776833351572e-27 },
        f64x2 { hi: -1.074122711280688e-10, lo: 4.648456344809929e-27 },
        f64x2 { hi: 4.795897620864847e-12, lo: 2.44769323702631e-28 },
        f64x2 { hi: 8.751221996380555e-13, lo: 3.070096873927715e-29 },
        f64x2 { hi: -2.1791367308069344e-14, lo: -1.1966032645603263e-30 },
        f64x2 { hi: -3.7357787089654046e-15, lo: 3.5273471362548918e-31 },
        f64x2 { hi: 2.19627024729479e-17, lo: -3.448285936206246e-34 },
        f64x2 { hi: 8.765871685007371e-18, lo: -6.275414546343881e-35 },
        f64x2 { hi: 5.724639270023456e-20, lo: -2.094006581883843e-36 },
        f64x2 { hi: -1.2161694920959997e-20, lo: 4.081072812178618e-38 },
        f64x2 { hi: -1.8684368042043668e-22, lo: -4.9606576573073464e-39 },
        f64x2 { hi: 1.0550891129477819e-23, lo: -5.5608324937398845e-40 },
        f64x2 { hi: 2.4861441687249707e-25, lo: -1.5010507446809888e-41 },
        f64x2 { hi: -5.823304748523291e-27, lo: 2.2079140848063186e-43 },
        f64x2 { hi: -2.0165183921069685e-28, lo: 6.992147028889838e-45 },
        f64x2 { hi: 1.896632671475821e-30, lo: 3.0052993839413775e-47 },
        f64x2 { hi: 1.1225095494047078e-31, lo: -8.682967183703833e-50 },
        f64x2 { hi: -1.6755789889608834e-34, lo: -7.828373160727297e-51 },
        f64x2 { hi: -4.543848034838703e-35, lo: 2.568853695440327e-51 },
    ],
    [
        f64x2 { hi: 1.228558508809108e-6, lo: -5.576438745979376e-23 },
        f64x2 { hi: -1.1940986396077243e-6, lo: 1.831606869872304e-23 },
        f64x2 { hi: -6.099999653919517e-8, lo: -4.5978238528928286e-24 },
        f64x2 { hi: -8.844063913885954e-9, lo: -2.8800816743411943e-25 },
        f64x2 { hi: 3.169816317194402e-8, lo: -1.397516653454564e-24 },
        f64x2 { hi: -1.4200472095883398e-8, lo: 4.224647971662184e-25 },
        f64x2 { hi: 3.161410591547148e-9, lo: -3.6795963344618467e-26 },
        f64x2 { hi: -2.443631526211608e-10, lo: -8.541055686027282e-27 },
        f64x2 { hi: -4.3226312365634374e-11, lo: -3.0356619165097975e-27 },
        f64x2 { hi: 9.017681907739495e-12, lo: 4.05861086639271e-28 },
        f64x2 { hi: 1.469890792000892e-13, lo: 1.0138544653952742e-29 },
        f64x2 { hi: -8.703305382470976e-14, lo: -4.094222013857739e-30 },
        f64x2 { hi: -8.379770803373182e-16, lo: -4.582637078399071e-32 },
        f64x2 { hi: 3.8874550686659373e-16, lo: 5.439052530648429e-34 },
        f64x2 { hi: 6.2406850724701105e-18, lo: 3.3302865220301215e-34 },
        f64x2 { hi: -9.22917087555887e-19, lo: 9.058575890854599e-35 },
        f64x2 { hi: -2.1592426398497923e-20, lo: 1.4454286166997812e-36 },
        f64x2 { hi: 1.264734812795467e-21, lo: 8.35237601026292e-38 },
        f64x2 { hi: 3.9909605998675126e-23, lo: -1.656511010571545e-39 },
        f64x2 { hi: -1.0363603911456208e-24, lo: 1.500882370275534e-41 },
        f64x2 { hi: -4.536116477268368e-26, lo: 6.409785275266073e-43 },
        f64x2 { hi: 4.764552198290305e-28, lo: 5.942039063871006e-45 },
        f64x2 { hi: 3.4582801513416345e-29, lo: -5.061282355925837e-48 },
        f64x2 { hi: -5.49948064974206e-32, lo: 9.14054126217011e-49 },
        f64x2 { hi: -1.8682482317645924e-32, lo: -6.135478418430844e-49 },
        f64x2 { hi: -8.445695356803679e-35, lo: 3.6446216632729875e-51 },
    ],
    [
        f64x2 { hi: -4.033411142597758e-6, lo: -1.212299158198647e-22 },
        f64x2 { hi: 2.0252281974869307e-6, lo: -6.512216921258346e-23 },
        f64x2 { hi: -6.11323732627802e-7, lo: 3.93331997773018e-23 },
        f64x2 { hi: 1.6899332657723014e-7, lo: 4.7337141697234086e-24 },
        f64x2 { hi: -3.8677374321150255e-8, lo: 3.4282078177559163e-25 },
        f64x2 { hi: 6.259906358926756e-9, lo: -3.966410796701389e-25 },
        f64x2 { hi: -3.6284669051529545e-10, lo: -1.1542800512367311e-26 },
        f64x2 { hi: -1.0805905023264905e-10, lo: -5.874360569251213e-28 },
        f64x2 { hi: 2.7038403322375307e-11, lo: 3.5850850920497545e-28 },
        f64x2 { hi: -1.225126787326331e-12, lo: 4.0447021056162443e-29 },
        f64x2 { hi: -2.7853879787768917e-13, lo: 1.387066711059111e-29 },
        f64x2 { hi: 2.21554370252174e-14, lo: 7.269602508830118e-31 },
        f64x2 { hi: 1.639404787914064e-15, lo: 2.931855551765088e-32 },
        f64x2 { hi: -1.1419338198194077e-16, lo: -6.817231154511536e-33 },
        f64x2 { hi: -6.4772080448473376e-18, lo: 1.9161470679487593e-34 },
        f64x2 { hi: 2.7686184008359666e-19, lo: 5.232466955064909e-36 },
        f64x2 { hi: 1.627436743160784e-20, lo: 1.2563794274341141e-36 },
        f64x2 { hi: -3.516446585292969e-22, lo: -2.527312545498901e-39 },
        f64x2 { hi: -2.6024300244972582e-23, lo: -1.4106867967454613e-39 },
        f64x2 { hi: 2.128054100563576e-25, lo: -4.288959557426396e-42 },
        f64x2 { hi: 2.7577285854639536e-26, lo: 1.15486767087008e-42 },
        f64x2 { hi: 1.496365963575927e-29, lo: 7.668894931563154e-46 },
        f64x2 { hi: -2.0211374502493027e-29, lo: 1.0510394975587957e-45 },
        f64x2 { hi: -1.3858418192488383e-31, lo: 8.153692638400925e-48 },
        f64x2 { hi: 1.0582814365469178e-32, lo: -3.1606671552556818e-49 },
        f64x2 { hi: 1.2971678942901192e-34, lo: 1.0091990020765233e-50 },
    ],
    [
        f64x2 { hi: 6.981157928224481e-8, lo: 3.9534655123334875e-24 },
        f64x2 { hi: 5.187602099781909e-8, lo: -3.2122050375061077e-24 },
        f64x2 { hi: -1.5025689400416704e-7, lo: 2.8399142863033162e-24 },
        f64x2 { hi: 5.385175415429129e-8, lo: 1.818414329268543e-24 },
        f64x2 { hi: -1.2009470947212667e-8, lo: 1.2315971496464968e-25 },
        f64x2 { hi: 1.8441416112134065e-9, lo: 1.538905835243495e-26 },
        f64x2 { hi: -6.051285922581879e-11, lo: -6.349866829419425e-27 },
        f64x2 { hi: -5.891392764479414e-11, lo: 1.4113169961978372e-27 },
        f64x2 { hi: 1.6515772641435116e-11, lo: 1.13572171248774e-27 },
        f64x2 { hi: -1.6489918275452742e-12, lo: -5.884341090021106e-29 },
        f64x2 { hi: -8.450007409241396e-14, lo: -5.618297160829461e-30 },
        f64x2 { hi: 3.023518017772655e-14, lo: 1.9392068642090218e-30 },
        f64x2 { hi: -6.17920112377458e-16, lo: -2.1017765124498416e-32 },
        f64x2 { hi: -2.1506480207808527e-16, lo: 2.5313177728609743e-33 },
        f64x2 { hi: 5.236058416945074e-18, lo: -2.4650203678575926e-34 },
        f64x2 { hi: 8.702944990758898e-19, lo: 6.90303316625067e-35 },
        f64x2 { hi: -1.2721127494561934e-20, lo: -9.917051222376512e-38 },
        f64x2 { hi: -2.1508806771393662e-21, lo: -8.245061847405083e-38 },
        f64x2 { hi: 9.442732617564068e-24, lo: 8.892996102795601e-41 },
        f64x2 { hi: 3.386918515792334e-24, lo: -6.676904911450353e-41 },
        f64x2 { hi: 1.2242358659889073e-26, lo: 6.554124299399133e-43 },
        f64x2 { hi: -3.5475860657436264e-27, lo: 1.7424162213646817e-43 },
        f64x2 { hi: -3.5154061138545925e-29, lo: -1.0967737227143288e-46 },
        f64x2 { hi: 2.5624450392466353e-30, lo: 1.3152414416639394e-46 },
        f64x2 { hi: 3.982565080027314e-32, lo: 2.1832998372564212e-48 },
        f64x2 { hi: -1.3043413384906499e-33, lo: 8.712760661195012e-51 },
    ],
];

pub(crate) const EXP_POLY_EXP_F64: [Complex<f64>; 18] = [
    Complex { re: 1.0, im: 1.0812144107799266e-23 },
    Complex { re: -4.479610786345225e-21, im: 1.0 },
    Complex { re: -0.5, im: 3.027415987199093e-19 },
    Complex { re: -8.413828845781633e-18, im: -0.16666666666666669 },
    Complex { re: 0.04166666666666679, im: 1.1656663448809618e-16 },
    Complex { re: -1.0603511152022404e-15, im: 0.008333333333332324 },
    Complex { re: -0.0013888888888827402, im: 5.789264486508273e-15 },
    Complex { re: -2.491995923872859e-14, im: -0.00019841269843586228 },
    Complex { re: 2.4801587374768556e-5, im: 6.704175576866034e-14 },
    Complex { re: -1.594987515102099e-13, im: 2.7557317787217356e-6 },
    Complex { re: -2.755729303110001e-7, im: 2.3127460502103687e-13 },
    Complex { re: -3.2663668749921504e-13, im: -2.5052389834713885e-8 },
    Complex { re: 2.087985316554709e-9, im: 2.5867211760028217e-13 },
    Complex { re: -2.2167241850689593e-13, im: 1.6041263496425594e-10 },
    Complex { re: -1.1352710114429515e-11, im: 8.943908448871146e-14 },
    Complex { re: -4.542339711641447e-14, im: -7.962911435347713e-13 },
    Complex { re: 5.979573239083729e-14, im: 7.185782517642856e-15 },
    Complex { re: -1.970149077208406e-15, im: 1.9701490772084063e-15 },
];
