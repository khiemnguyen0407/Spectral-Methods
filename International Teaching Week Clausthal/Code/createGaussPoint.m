function [gW, gP] = createGaussPoint(numGauss)
%CREATEGAUSSPOINT  Gauss points and the associated weights.

switch numGauss
    case 1
        gP = 0;
        gW = 2;
    case 2
        gP = [-1/sqrt(3)  1/sqrt(3)];
        gW = [1  1];
    case 3
        gP = [-sqrt(3/5)   0   sqrt(3/5)];
        gW = [5/9   8/9   5/9];
    case 4
        a = [-0.86113631159405257522394648889281, ...
            -0.33998104358485626480266575910324];
        b = [0.34785484513745385737306394922200, ...
            0.65214515486254614262693605077800];
        gP = horzcat(a,fliplr(-a));
        gW = horzcat(b,fliplr(b));
    case 5
        a = [-0.90617984593866399279762687829939, ...
            -0.53846931010568309103631442070021];
        b = [0.23692688505618908751426404071992, ...
            0.47862867049936646804129151483564];
        mp = 0;
        mw = 0.56888888888888888888888888888889;
        gP = horzcat(a,mp,fliplr(-a));
        gW = horzcat(b,mw,fliplr(b));
    case 6
        a = [-0.93246951420315202781230155449399, ...
            -0.66120938646626451366139959501991, ...
            -0.23861918608319690863050172168071];
        b = [0.17132449237917034504029614217273, ...
            0.36076157304813860756983351383772, ...
            0.46791393457269104738987034398955];
        gP = horzcat(a,fliplr(-a));
        gW = horzcat(b,fliplr(b));
    case 7
        a = [-0.94910791234275852452618968404785, ...
            -0.74153118559939443986386477328079, ...
            -0.40584515137739716690660641207696];
        b = [0.12948496616886969327061143267908, ...
            0.27970539148927666790146777142378, ...
            0.38183005050511894495036977548898];
        mp = 0;
        mw = 0.41795918367346938775510204081633;
        gP = horzcat(a,mp,fliplr(-a));
        gW = horzcat(b,mw,fliplr(b));
    case 8
        a = [-0.96028985649753623168356086856947, ...
            -0.79666647741362673959155393647583, ...
            -0.52553240991632898581773904918925, ...
            -0.18343464249564980493947614236018];
        b = [0.10122853629037625915253135430996, ...
            0.22238103445337447054435599442624, ...
            0.31370664587788728733796220198660, ...
            0.36268378337836198296515044927720];
        gP = horzcat(a,fliplr(-a));
        gW = horzcat(b,fliplr(b));
    case 9
        a = [-0.96816023950762608983557620290367, ...
            -0.83603110732663579429942978806973, ...
            -0.61337143270059039730870203934147, ...
            -0.32425342340380892903853801464334];
        b = [0.081274388361574411971892158110524, ...
            0.18064816069485740405847203124291, ...
            0.26061069640293546231874286941863, ...
            0.31234707704000284006863040658444];
        mp = 0;
        mw = 0.33023935500125976316452506928697;
        gP = horzcat(a,mp,fliplr(-a));
        gW = horzcat(b,mw,fliplr(b));
    case 10
        a = [-0.97390652851717172007796401208445, ...
            -0.86506336668898451073209668842349, ...
            -0.67940956829902440623432736511487, ...
            -0.43339539412924719079926594316578, ...
            -0.14887433898163121088482600112972];
        b = [0.066671344308688137593568809893332, ...
            0.14945134915058059314577633965770, ...
            0.21908636251598204399553493422816, ...
            0.26926671930999635509122692156947, ...
            0.29552422471475287017389299465134];
        gP = horzcat(a,fliplr(-a));
        gW = horzcat(b,fliplr(b));
    case 11
        a = [-0.97822865814605699280393800112286, ...
            -0.88706259976809529907515776930393, ...
            -0.73015200557404932409341625203115, ...
            -0.51909612920681181592572566945861, ...
            -0.26954315595234497233153198540086];
        b = [0.055668567116173666482753720442549, ...
            0.12558036946490462463469429922394, ...
            0.18629021092773425142609764143166, ...
            0.23319376459199047991852370484318, ...
            0.26280454451024666218068886989051];
        mp = 0;
        mw = 0.27292508677790063071448352833634;
        gP = horzcat(a,mp,fliplr(-a));
        gW = horzcat(b,mw,fliplr(b));
    case 12
        a = [-0.98156063424671925069054909014928, ...
            -0.90411725637047485667846586611910, ...
            -0.76990267419430468703689383321282, ...
            -0.58731795428661744729670241894053, ...
            -0.36783149899818019375269153664372, ...
            -0.12523340851146891547244136946385];
        b = [0.047175336386511827194615961485017, ...
            0.10693932599531843096025471819400, ...
            0.16007832854334622633465252954336, ...
            0.20316742672306592174906445580980, ...
            0.23349253653835480876084989892488, ...
            0.24914704581340278500056243604295];
        gP = horzcat(a,fliplr(-a));
        gW = horzcat(b,fliplr(b));
    case 13
        a = [-0.98418305471858814947282944880711, ...
            -0.91759839922297796520654783650072, ...
            -0.80157809073330991279420648958286, ...
            -0.64234933944034022064398460699552, ...
            -0.44849275103644685287791285212764, ...
            -0.23045831595513479406552812109799];
        b = [0.040484004765315879520021592200986, ...
            0.092121499837728447914421775953797, ...
            0.13887351021978723846360177686887, ...
            0.17814598076194573828004669199610, ...
            0.20781604753688850231252321930605, ...
            0.22628318026289723841209018603978];
        mp = 0;
        mw = 0.23255155323087391019458951526884;
        gP = horzcat(a,mp,fliplr(-a));
        gW = horzcat(b,mw,fliplr(b));
    case 14
        a = [-0.98628380869681233884159726670405, ...
            -0.92843488366357351733639113937787, ...
            -0.82720131506976499318979474265039, ...
            -0.68729290481168547014801980301933, ...
            -0.51524863635815409196529071855119, ...
            -0.31911236892788976043567182416848, ...
            -0.10805494870734366206624465021983];
        b = [0.035119460331751863031832876138192, ...
            0.080158087159760209805633277062854, ...
            0.12151857068790318468941480907248, ...
            0.15720316715819353456960193862384, ...
            0.18553839747793781374171659012516, ...
            0.20519846372129560396592406566122, ...
            0.21526385346315779019587644331626];
        gP = horzcat(a,fliplr(-a));
        gW = horzcat(b,fliplr(b));
    case 15
        a = [-0.98799251802048542848956571858661, ...
            -0.93727339240070590430775894771021, ...
            -0.84820658341042721620064832077422, ...
            -0.72441773136017004741618605461394, ...
            -0.57097217260853884753722673725391, ...
            -0.39415134707756336989720737098105, ...
            -0.20119409399743452230062830339460];
        b = [0.030753241996117268354628393577204, ...
            0.070366047488108124709267416450667, ...
            0.10715922046717193501186954668587, ...
            0.13957067792615431444780479451103, ...
            0.16626920581699393355320086048121, ...
            0.18616100001556221102680056186642, ...
            0.19843148532711157645611832644384];
        mp = 0;
        mw = 0.20257824192556127288062019996752;
        gP = horzcat(a,mp,fliplr(-a));
        gW = horzcat(b,mw,fliplr(b));
    case 16
        a = [-0.98940093499164993259615417345033, ...
            -0.94457502307323257607798841553461, ...
            -0.86563120238783174388046789771239, ...
            -0.75540440835500303389510119484744, ...
            -0.61787624440264374844667176404879, ...
            -0.45801677765722738634241944298358, ...
            -0.28160355077925891323046050146050, ...
            -0.095012509837637440185319335424958];
        b = [0.027152459411754094851780572456018, ...
            0.062253523938647892862843836994378, ...
            0.095158511682492784809925107602246, ...
            0.12462897125553387205247628219202, ...
            0.14959598881657673208150173054748, ...
            0.16915651939500253818931207903036, ...
            0.18260341504492358886676366796922, ...
            0.18945061045506849628539672320828];
        gP = horzcat(a,fliplr(-a));
        gW = horzcat(b,fliplr(b));
    case 17
        a = [-0.99057547531441733567543401994067, ...
            -0.95067552176876776122271695789580, ...
            -0.88023915372698590212295569448816, ...
            -0.78151400389680140692523005552048, ...
            -0.65767115921669076585030221664300, ...
            -0.51269053708647696788624656862955, ...
            -0.35123176345387631529718551709535, ...
            -0.17848418149584785585067749365407];
        b = [0.024148302868547931960110026287565, ...
            0.055459529373987201129440165358245, ...
            0.085036148317179180883535370191062, ...
            0.11188384719340397109478838562636, ...
            0.13513636846852547328631998170235, ...
            0.15404576107681028808143159480196, ...
            0.16800410215645004450997066378832, ...
            0.17656270536699264632527099011320];
        mp = 0;
        mw = 0.17944647035620652545826564426189;
        gP = horzcat(a,mp,fliplr(-a));
        gW = horzcat(b,mw,fliplr(b));
    case 18
        a = [-0.99156516842093094673001600470615, ...
            -0.95582394957139775518119589292978, ...
            -0.89260246649755573920606059112715, ...
            -0.80370495897252311568241745501459, ...
            -0.69168704306035320787489108128885, ...
            -0.55977083107394753460787154852533, ...
            -0.41175116146284264603593179383305, ...
            -0.25188622569150550958897285487791, ...
            -0.084775013041735301242261852935784];
        b = [0.021616013526483310313342710266452, ...
            0.049714548894969796453334946202639, ...
            0.076425730254889056529129677616637, ...
            0.10094204410628716556281398492483, ...
            0.12255520671147846018451912680020, ...
            0.14064291467065065120473130375195, ...
            0.15468467512626524492541800383637, ...
            0.16427648374583272298605377646593, ...
            0.16914238296314359184065647013499];
        gP = horzcat(a,fliplr(-a));
        gW = horzcat(b,fliplr(b));
    case 19
        a = [-0.99240684384358440318901767025326, ...
            -0.96020815213483003085277884068765, ...
            -0.90315590361481790164266092853231, ...
            -0.82271465653714282497892248671271, ...
            -0.72096617733522937861709586082378, ...
            -0.60054530466168102346963816494624, ...
            -0.46457074137596094571726714810410, ...
            -0.31656409996362983199011732884984, ...
            -0.16035864564022537586809611574074];
        b = [0.019461788229726477036312041464438, ...
            0.044814226765699600332838157401994, ...
            0.069044542737641226580708258006013, ...
            0.091490021622449999464462094123840, ...
            0.11156664554733399471602390168177, ...
            0.12875396253933622767551578485688, ...
            0.14260670217360661177574610944190, ...
            0.15276604206585966677885540089766, ...
            0.15896884339395434764995643946505];
        mp = 0;
        mw = 0.16105444984878369597916362532092;
        gP = horzcat(a,mp,fliplr(-a));
        gW = horzcat(b,mw,fliplr(b));
    case 20
        a = [-0.99312859918509492478612238847132, ...
            -0.96397192727791379126766613119728, ...
            -0.91223442825132590586775244120330, ...
            -0.83911697182221882339452906170152, ...
            -0.74633190646015079261430507035564, ...
            -0.63605368072651502545283669622629, ...
            -0.51086700195082709800436405095525, ...
            -0.37370608871541956067254817702493, ...
            -0.22778585114164507808049619536857, ...
            -0.076526521133497333754640409398838];
        b = [0.017614007139152118311861962351853, ...
            0.040601429800386941331039952274932, ...
            0.062672048334109063569506535187042, ...
            0.083276741576704748724758143222046, ...
            0.10193011981724043503675013548035, ...
            0.11819453196151841731237737771138, ...
            0.13168863844917662689849449974816, ...
            0.14209610931838205132929832506716, ...
            0.14917298647260374678782873700197, ...
            0.15275338713072585069808433195510];
        gP = horzcat(a,fliplr(-a));
        gW = horzcat(b,fliplr(b));
    case 21
        a = [-0.99375217062038950026024203593794, ...
            -0.96722683856630629431662221490770, ...
            -0.92009933415040082879018713371497, ...
            -0.85336336458331728364725063858757, ...
            -0.76843996347567790861587785130623, ...
            -0.66713880419741231930596666999034, ...
            -0.55161883588721980705901879672431, ...
            -0.42434212020743878357366888854379, ...
            -0.28802131680240109660079251606460, ...
            -0.14556185416089509093703098233869];
        b = [0.016017228257774333324224616858471, ...
            0.036953789770852493799950668299330, ...
            0.057134425426857208283635826472448, ...
            0.076100113628379302017051653300183, ...
            0.093444423456033861553289741113932, ...
            0.10879729916714837766347457807011, ...
            0.12183141605372853419536717712573, ...
            0.13226893863333746178105257449678, ...
            0.13988739479107315472213342386758, ...
            0.14452440398997005906382716655375];
        mp = 0;
        mw = 0.14608113364969042719198514768337;
        gP = horzcat(a,mp,fliplr(-a));
        gW = horzcat(b,mw,fliplr(b));
    case 22
        a = [-0.99429458548239929207303142116130, ...
            -0.97006049783542872712395098676527, ...
            -0.92695677218717400052069293925905, ...
            -0.86581257772030013653642563701938, ...
            -0.78781680597920816200427795540835, ...
            -0.69448726318668278005068983576226, ...
            -0.58764040350691159295887692763865, ...
            -0.46935583798675702640633071096641, ...
            -0.34193582089208422515814742042738, ...
            -0.20786042668822128547884653391955, ...
            -0.069739273319722221213841796118628];
        b = [0.014627995298272200684991098047185, ...
            0.033774901584814154793302246865913, ...
            0.052293335152683285940312051273211, ...
            0.069796468424520488094961418930218, ...
            0.085941606217067727414443681372703, ...
            0.10041414444288096493207883783054, ...
            0.11293229608053921839340060742178, ...
            0.12325237681051242428556098615481, ...
            0.13117350478706237073296499253031, ...
            0.13654149834601517135257383123152, ...
            0.13925187285563199337541024834181];
        gP = horzcat(a,fliplr(-a));
        gW = horzcat(b,fliplr(b));
    case 23
        a = [-0.99476933499755212352392571544557, ...
            -0.97254247121811523195602407682078, ...
            -0.93297108682601610234919698903842, ...
            -0.87675235827044166737815688593415, ...
            -0.80488840161883989215111840699678, ...
            -0.71866136313195019446162448374862, ...
            -0.61960987576364615638509731164960, ...
            -0.50950147784600754968979304786685, ...
            -0.39030103803029083142148887288061, ...
            -0.26413568097034493053386953828331, ...
            -0.13325682429846611093174268224177];
        b = [0.013411859487141772081309493458615, ...
            0.030988005856979444310694219641885, ...
            0.048037671731084668571641071632034, ...
            0.064232421408525852127169615158911, ...
            0.079281411776718954922892524742043, ...
            0.092915766060035147477018617369765, ...
            0.10489209146454141007408618501474, ...
            0.11499664022241136494164351293396, ...
            0.12304908430672953046757840067201, ...
            0.12890572218808214997859533939979, ...
            0.13246203940469661737164246470332];
        mp = 0;
        mw = 0.13365457218610617535145711054584;
        gP = horzcat(a,mp,fliplr(-a));
        gW = horzcat(b,mw,fliplr(b));
    case 24
        a = [-0.99518721999702136017999740970074, ...
            -0.97472855597130949819839199300817, ...
            -0.93827455200273275852364900170872, ...
            -0.88641552700440103421315434198220, ...
            -0.82000198597390292195394987266975, ...
            -0.74012419157855436424382810309998, ...
            -0.64809365193697556925249578691075, ...
            -0.54542147138883953565837561721837, ...
            -0.43379350762604513848708423191335, ...
            -0.31504267969616337438679329131981, ...
            -0.19111886747361630915863982075707, ...
            -0.064056892862605626085043082624745];
        b = [0.012341229799987199546805667070037, ...
            0.028531388628933663181307815951878, ...
            0.044277438817419806168602748211338, ...
            0.059298584915436780746367758500109, ...
            0.073346481411080305734033615253117, ...
            0.086190161531953275917185202983743, ...
            0.097618652104113888269880664464247, ...
            0.10744427011596563478257734244661, ...
            0.11550566805372560135334448390678, ...
            0.12167047292780339120446315347626, ...
            0.12583745634682829612137538251118, ...
            0.12793819534675215697405616522470];
        gP = horzcat(a,fliplr(-a));
        gW = horzcat(b,fliplr(b));
    case 25
    case 26
    case 27
    case 28
    case 29
    case 30
    case 31
    case 32
        a = [-0.99726386184948156354498112866504, ...
            -0.98561151154526833540017504463090, ...
            -0.96476225558750643077381192811827, ...
            -0.93490607593773968917091913483541, ...
            -0.89632115576605212396530724371921, ...
            -0.84936761373256997013369300496774, ...
            -0.79448379596794240696309729897043, ...
            -0.73218211874028968038742666509127, ...
            -0.66304426693021520097511516866324, ...
            -0.58771575724076232904074547640183, ...
            -0.50689990893222939002374747437782, ...
            -0.42135127613063534536411943617243, ...
            -0.33186860228212764977991680573019, ...
            -0.23928736225213707454460320916550, ...
            -0.14447196158279649348518637359881, ...
            -0.048307665687738316234812570440502];
        b = [0.0070186100094700966004070637388532, ...
            0.016274394730905670605170562206387, ...
            0.025392065309262059455752589789224, ...
            0.034273862913021433102687732252373, ...
            0.042835898022226680656878646606126, ...
            0.050998059262376176196163244689522, ...
            0.058684093478535547145283637300171, ...
            0.065822222776361846837650063706939, ...
            0.072345794108848506225399356478488, ...
            0.078193895787070306471740918828307, ...
            0.083311924226946755222199074604349, ...
            0.087652093004403811142771462751802, ...
            0.091173878695763884712868577111637, ...
            0.093844399080804565639180237668117, ...
            0.095638720079274859419082002204131, ...
            0.096540088514727800566764830063576];
        gP = horzcat(a,fliplr(-a));
        gW = horzcat(b,fliplr(b));
    case 33
    case 34
    case 35
    case 36
    case 37
    case 38
    case 39
    case 40
    case 41
    case 42
    case 43
    case 44
    case 45
    case 46
    case 47
    case 48
    case 49
    case 50
    case 51
    case 52
    case 53
    case 54
    case 55
    case 56
    case 57
    case 58
    case 59
    case 60
    case 61
    case 62
    case 63
    case 64
        a = [-0.99930504173577213945690562434564, ...
            -0.99634011677195527934692450067640, ...
            -0.99101337147674432073938238344330, ...
            -0.98333625388462595693129930215683, ...
            -0.97332682778991096374185350735227, ...
            -0.96100879965205371891861412189716, ...
            -0.94641137485840281606248149134726, ...
            -0.92956917213193957582149015455923, ...
            -0.91052213707850280575638066800833, ...
            -0.88931544599511410585340403827285, ...
            -0.86599939815409281976078338507016, ...
            -0.84062929625258036275169154469587, ...
            -0.81326531512279755974192333808630, ...
            -0.78397235894334140761022052521377, ...
            -0.75281990726053189661186377488569, ...
            -0.71988185017161082684894021783195, ...
            -0.68523631305423324256355837103138, ...
            -0.64896547125465733985776123199340, ...
            -0.61115535517239325024885297101855, ...
            -0.57189564620263403428387811665919, ...
            -0.53127946401989454565801390354446, ...
            -0.48940314570705295747852630702192, ...
            -0.44636601725346408798494771475892, ...
            -0.40227015796399160369576677126016, ...
            -0.35722015833766811595044261504620, ...
            -0.31132287199021095615751269856016, ...
            -0.26468716220876741637396417251002, ...
            -0.21742364374000708414964874898882, ...
            -0.16964442042399281803731362974827, ...
            -0.12146281929612055447037646349225, ...
            -0.072993121787799039449542941940337, ...
            -0.024350292663424432508955842853716];
        b = [0.0017832807216964329472960791449719, ...
            0.0041470332605624676352875357285514, ...
            0.0065044579689783628561173603999813, ...
            0.0088467598263639477230309146597306, ...
            0.011168139460131128818590493019208, ...
            0.013463047896718642598060766685956, ...
            0.015726030476024719321965995297540, ...
            0.017951715775697343085045302001119, ...
            0.020134823153530209372340316728544, ...
            0.022270173808383254159298330384155, ...
            0.024352702568710873338177550409069, ...
            0.026377469715054658671691792625225, ...
            0.028339672614259483227511305200237, ...
            0.030234657072402478867974059819549, ...
            0.032057928354851553585467504347899, ...
            0.033805161837141609391565482110725, ...
            0.035472213256882383810693146715246, ...
            0.037055128540240046040415101809583, ...
            0.038550153178615629128962496946809, ...
            0.039953741132720341386656926128336, ...
            0.041262563242623528610156297473638, ...
            0.042473515123653589007339767908817, ...
            0.043583724529323453376827860973737, ...
            0.044590558163756563060134710030945, ...
            0.045491627927418144479770996971269, ...
            0.046284796581314417295953249232261, ...
            0.046968182816210017325326285754581, ...
            0.047540165714830308662282206944223, ...
            0.047999388596458307728126179871346, ...
            0.048344762234802957169769527158018, ...
            0.048575467441503426934799066783978, ...
            0.048690957009139720383365390734750];
        gP = horzcat(a,fliplr(-a));
        gW = horzcat(b,fliplr(b));
    otherwise
        gP = 0;
        gW = 2;
end