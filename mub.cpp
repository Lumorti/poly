#include "poly.h"
#include <time.h>

// Standard cpp entry point 
int main(int argc, char ** argv) {

	// Sometimes we want this, sometimes we don't
	std::srand(std::time(NULL));

	// Get the problem from the args
	int d = 2;
	int n = 2;
	if (argc > 1) {
		d = std::stoi(argv[1]);
	}
	if (argc > 2) {
		n = std::stoi(argv[2]);
	}

	// Useful quantities
	int numVarsNonConj = n*d*d;
	int numVars = 2*numVarsNonConj;
	int conjDelta = numVarsNonConj;
	double rt2 = 1.0/std::sqrt(2.0);

	// The list of equations to fill
	std::vector<Polynomial<double>> eqns;

	// Generate equations
	std::cout << "Generating equations..." << std::endl;
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			for (int k=0; k<d; k++) {
				for (int l=0; l<d; l++) {

					// (a+ib)*(c+id)
					Polynomial<std::complex<double>> eqn(numVars);
					for (int m=0; m<d; m++) {
						int var1 = i*d*d + k*d + m;
						int var2 = j*d*d + l*d + m;
						int var3 = var1 + conjDelta;
						int var4 = var2 + conjDelta;
						eqn.addTerm(1, {var1, var2});
						eqn.addTerm(1, {var3, var4});
						eqn.addTerm(1i, {var1, var4});
						eqn.addTerm(-1i, {var2, var3});
					}

					// For the normalisation equations
					if (i == j && k == l) {
						eqn.addTerm(-1.0, {});
					}

					// For the MUB-ness equations
					if (i != j) {
						eqn = eqn.conjugate() * eqn;
						eqn.addTerm(-1.0/d, {});
	 				}

					// Both the real and imag parts should be 0
					eqns.push_back(Polynomial<double>(eqn.real()));
					eqns.push_back(Polynomial<double>(eqn.imag()));

				}
			}
		}
	}

	// For replacing variables with known ideal values
	std::vector<Polynomial<double>> newEqns;
	std::vector<int> indsToReplace = {};
	std::vector<double> valsToReplaceReal = {};
	std::vector<double> valsToReplaceImag = {};
	std::vector<double> toAdd = {};

	// For d == 2
	if (d == 2) {
		if (n >= 2) {
			toAdd = {1, 0, 0, 1, 0, 0, 0, 0};
			valsToReplaceReal.insert(valsToReplaceReal.end(), toAdd.begin(), toAdd.begin()+d*d);
			valsToReplaceImag.insert(valsToReplaceImag.end(), toAdd.begin()+d*d, toAdd.end());
		} 
		if (n >= 3) {
			//toAdd = {0.32195961348437352, -0.36877176426602609, 0.28247707126056093, -0.70033770335114875, -0.62955699833396306, -0.60333037184645466, 0.64823358649198259, -0.097606780602686141};
			//valsToReplaceReal.insert(valsToReplaceReal.end(), toAdd.begin(), toAdd.begin()+d*d);
			//valsToReplaceImag.insert(valsToReplaceImag.end(), toAdd.begin()+d*d, toAdd.end());
		} 
		if (n >= 4) {
			//toAdd = {-0.018102279444082173, -0.38456996779632113, 0.70685839058825606, -0.61261003259874924, 0.7068750491645367, 0.59338519782444576, -0.01873958206243807, -0.35314161801162253};
			//valsToReplaceReal.insert(valsToReplaceReal.end(), toAdd.begin(), toAdd.begin()+d*d);
			//valsToReplaceImag.insert(valsToReplaceImag.end(), toAdd.begin()+d*d, toAdd.end());
		}
	}

	// For d == 3
	if (d == 3) {
		if (n >= 2) {
			toAdd = {1, 0, 0, 0, 1, 0, 0, 0, 1,  0,0,0,0,0,0,0,0,0};
			valsToReplaceReal.insert(valsToReplaceReal.end(), toAdd.begin(), toAdd.begin()+d*d);
			valsToReplaceImag.insert(valsToReplaceImag.end(), toAdd.begin()+d*d, toAdd.end());
		} 
		if (n >= 3) {
			//toAdd = {0.38003647513790584, 0.54944290475413904, 0.56175009870872716, -0.53097781219693641, 0.57425807822404973, -0.27641787629598696, -0.55668556729245144, -0.16226138377826449, 0.57002016182175375, 0.43463273300760596, 0.17732972525956039, 0.13330481838767158, 0.22670670859154296, 0.059674045994163538, -0.50687918151935973, 0.15308340614683558, -0.55407996643952395, 0.091707902369081565};
			//valsToReplaceReal.insert(valsToReplaceReal.end(), toAdd.begin(), toAdd.begin()+d*d);
			//valsToReplaceImag.insert(valsToReplaceImag.end(), toAdd.begin()+d*d, toAdd.end());
		} 
		if (n >= 4) {
			//toAdd = {-0.42944516596525989, 0.56280358350677828, 0.5508507084147074, 0.33542576379070005, 0.52932788713914003, -0.43570786463396993, 0.53368684005965727, -0.21140477642105621, 0.56242000608959197, 0.38588887634933783, -0.12878450558033511, -0.17290694570063103, 0.46991793689607042, 0.23053269628828529, 0.37880341643143572, 0.22025368599321965, 0.53725354162382766, -0.13044949094764793};
			//valsToReplaceReal.insert(valsToReplaceReal.end(), toAdd.begin(), toAdd.begin()+d*d);
			//valsToReplaceImag.insert(valsToReplaceImag.end(), toAdd.begin()+d*d, toAdd.end());
		}
		if (n >= 5) {
			//toAdd = {-0.56239763060056369, 0.55502663659237716, -0.38216528000855249, -0.24367163760443913, 0.060175507569309686, 0.43888755922401024, -0.27135268664392387, 0.50505520785677116, 0.48134085556057499, 0.13054591686873751, 0.15899296029218193, 0.43276209924083153, 0.52340948625702755, 0.57420577021127073, -0.37511467578255853, -0.50960874569864689, -0.27973658913525329, -0.31881709742397063};
			//valsToReplaceReal.insert(valsToReplaceReal.end(), toAdd.begin(), toAdd.begin()+d*d);
			//valsToReplaceImag.insert(valsToReplaceImag.end(), toAdd.begin()+d*d, toAdd.end());
		}

	}

	// For d == 4 TODO d4n4
	if (d == 4) {
		if (n >= 2) {
			toAdd = {1, 0, 0, 0,   0, 1, 0, 0,    0, 0, 1, 0,   0, 0, 0, 1,     0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0};
			valsToReplaceReal.insert(valsToReplaceReal.end(), toAdd.begin(), toAdd.begin()+d*d);
			valsToReplaceImag.insert(valsToReplaceImag.end(), toAdd.begin()+d*d, toAdd.end());
		} 
		if (n >= 3) {
			//toAdd = {-0.16455283827451397, 0.46357286733819097, -0.22693280927043485, 0.32748058852088419, 0.46904250273404152, 0.43936237341855411, -0.1995756280179663, 0.49984251822243414, -0.46428216110854148, 0.11354695463316343, -0.18734714317541029, 0.12055068772022119, -0.49977240818388696, -0.29837475395247037, -0.35672793988614476, 0.083041662723964335, -0.47214628616943571, 0.18734939678452348, 0.44553562636143679, 0.37783288330025488, 0.17320702323392326, -0.23866680666752341, -0.45844322832185264, -0.01251124978706007, -0.18558998340793553, 0.48693690530779227, -0.46357386723882943, -0.48525021077864716, 0.015036932594357228, -0.40121421853279643, -0.35035245228450973, 0.49305764417729936};
			//valsToReplaceReal.insert(valsToReplaceReal.end(), toAdd.begin(), toAdd.begin()+d*d);
			//valsToReplaceImag.insert(valsToReplaceImag.end(), toAdd.begin()+d*d, toAdd.end());
		} 
		if (n >= 4) {
			//toAdd = {-0.34517723490886731, -0.46438772102691367, 0.4162772474269561, 0.32570711905589406, 0.18120232076003334, -0.35866513057479699, 0.49244149336385368, 0.025861072682707788, 0.494562438949945, 0.4191251201424997, 0.46628325440335683, -0.12224038228314307, 0.44480612732193858, -0.30946480875201704, 0.13700172247646483, -0.49399035977648414, -0.11944188380213282, 0.093398188152192793, 0.4901642316868231, 0.038124753490016192, -0.4912404666488614, -0.40793175271257087, -0.29350182716588108, 0.43758751208913116, -0.47356518001013537, 0.36388142060891521, 0.49984458318995323, 0.048633001438790083, -0.16478459020551803, -0.046566943121783347, -0.29760988069796329, 0.36245486005184352, -0.36173564517197998, -0.18532146010790818, 0.27697157502415154, -0.37936113845948549, 0.46601042667587178, -0.34836665302966346, -0.086610477687320761, 0.49933075890637929, -0.073539078984496922, -0.27264286726687953, -0.1804991076052116, -0.48482707012900678, -0.2283582898570653, 0.39272322785663816, 0.48086435597073207, 0.077288571736255432, 0.48552408404413233, 0.4911993277905215, -0.098686502215763389, 0.49854438420108127, 0.093181559114212817, 0.28912226771678678, -0.40479214152785437, -0.24190322443735834, 0.16042450311951459, -0.34291443830275858, -0.012465642161714767, -0.4976292106027374, 0.4720657153092469, -0.49782679583287415, 0.40178148296303012, 0.34442194231495865};
			//valsToReplaceReal.insert(valsToReplaceReal.end(), toAdd.begin(), toAdd.begin()+d*d);
			//valsToReplaceImag.insert(valsToReplaceImag.end(), toAdd.begin()+d*d, toAdd.end());
		}
		if (n >= 5) {
			//toAdd = {};
			//valsToReplaceReal.insert(valsToReplaceReal.end(), toAdd.begin(), toAdd.begin()+d*d);
			//valsToReplaceImag.insert(valsToReplaceImag.end(), toAdd.begin()+d*d, toAdd.end());
		}
		if (n >= 6) {
			//toAdd = {};
			//valsToReplaceReal.insert(valsToReplaceReal.end(), toAdd.begin(), toAdd.begin()+d*d);
			//valsToReplaceImag.insert(valsToReplaceImag.end(), toAdd.begin()+d*d, toAdd.end());
		}

	}

	// For d == 6
	if (d == 6) {
		if (n >= 2) {
			toAdd = {1,0,0,0,0,0,  0,1,0,0,0,0,  0,0,1,0,0,0,  0,0,0,1,0,0,  0,0,0,0,1,0,  0,0,0,0,0,1,  0,0,0,0,0,0,  0,0,0,0,0,0,  0,0,0,0,0,0,  0,0,0,0,0,0,  0,0,0,0,0,0,  0,0,0,0,0,0};
			valsToReplaceReal.insert(valsToReplaceReal.end(), toAdd.begin(), toAdd.begin()+d*d);
			valsToReplaceImag.insert(valsToReplaceImag.end(), toAdd.begin()+d*d, toAdd.end());
		} 
		if (n >= 3) {
			//toAdd = {0.39274799312809788, -0.39939857127662759, -0.035489394817327248, -0.2236038463601667, 0.22458226851873636, -0.38281471569592884, 0.13789507884557656, -0.088350685252632871, -0.20000814398580366, -0.28602386887467462, -0.10486261802780397, 0.21632665008582555, -0.34670174199703568, -0.40808070071280084, -0.39224982906711919, 0.021926846846029274, -0.40463234235863454, 0.19553562430962476, -0.26589849162638407, 0.30560911273111013, -0.022881520292220368, -0.33248206111875978, -0.02780982680107999, -0.13746107809643435, 0.31627908528829241, -0.16593365337750912, -0.21179563451780672, 0.27733562009779761, -0.39984430492136214, 0.26011713572465167, 0.064631078477285345, 0.40704587822940808, 0.32954122726986107, -0.36813603123947197, -0.33565341476283028, 0.23682174652126373, -0.11142567433602517, -0.084542494052860795, -0.40670280429723793, -0.34156696336236414, 0.34092445610854255, -0.14184350473464194, -0.38425461242224623, -0.39857352126391976, 0.35589812776115876, -0.29130236087548012, -0.3945509571483139, -0.34622159176820938, 0.21555644211117439, -0.011698271349969171, -0.11316683716589625, -0.40765905805286534, 0.054215669931191257, 0.35837475321219875, 0.30978173030835676, 0.27068381508736267, -0.40760647455568605, -0.23690158067646699, -0.40730000006503986, -0.38441017618762324, -0.25813593576194455, 0.37300497604952332, -0.34901188767118158, 0.29958574760847501, -0.082409021578354752, -0.3146517841591891, -0.40309987970004185, 0.031309424356401834, -0.24097580397849563, -0.17647228533267995, 0.23238646283695744, 0.3325389511286645};
			//valsToReplaceReal.insert(valsToReplaceReal.end(), toAdd.begin(), toAdd.begin()+d*d);
			//valsToReplaceImag.insert(valsToReplaceImag.end(), toAdd.begin()+d*d, toAdd.end());
		} 
		if (n >= 4) {
			toAdd = {0.29851718316866716, 0.28351026641716154, 0.39612093911012164, 0.38211687082289814, -0.0080723504907833331, -0.3258541267504797, 0.35148772338881074, 0.059319210887180465, 0.20156314273815906, 0.35954495020892535, -0.37447007587921954, 0.25896472478617866, -0.101099390455981, -0.15397814126987686, 0.29031641737141151, 0.109309677402671, 0.051697151525550745, 0.33146175193542049, 0.4041893068116304, 0.29077259104258557, -0.37606574663065462, -0.32190109034918635, -0.35462024150520377, -0.34674620051749316, -0.2638231551995161, -0.3070288523982409, -0.079626296107071129, 0.12519288236359208, -0.40730938702360442, -0.26130355221258261, -0.03948113605557637, 0.36122854877079458, -0.40559264154155361, 0.36816080054390316, -0.089424334349845702, 0.10225558594726765, -0.38575287378430362, -0.40824599391663718, -0.36476774182965588, -0.40178848112675325, -0.39660369621897762, -0.31281652677877875, -0.40639245172501087, 0.39700589531422315, -0.24310060095085897, 0.1400885898022895, -0.073992071149250674, -0.11131594450696203, -0.020247716965854159, -0.18043358970077727, -0.0044212442786832565, 0.33487346020970471, 0.23804946850507724, -0.17469652307760267, -0.16070107388677207, -0.32281145244684312, 0.21582543735985832, -0.40790436192118384, 0.40061578081923666, -0.2629223124378382, 0.067636214174507575, 0.38963280760597579, 0.26480200485408351, -0.37714397854625048, 0.39522107313493554, -0.38768309960377123, -0.29237927098189587, -0.12177261721028194, 0.02125030245907094, 0.39443968179146716, 0.070767553768784769, -0.15816434062642942, 0.27848816431861206, -0.29374869760351224, 0.098765739122610205, 0.14371364902758488, -0.40816807589485288, 0.24593762042808823, -0.20765894501396817, 0.4039162555569773, 0.35501981281907269, 0.19337585599868845, 0.16259992721142072, -0.31560175911044153, 0.3955317946773883, -0.37809620100429386, -0.28702458741375575, 0.39334260928899045, 0.40496259585324046, -0.23832752871359944, -0.057423663528933036, -0.28656226885276459, -0.15887648797968557, 0.25109130596533624, 0.20226395261622698, -0.21548465011426571, 0.31155013087258387, -0.26907493091966489, 0.40040759549020588, -0.38857844024297922, 0.027667846941405318, -0.31366734680809394, 0.40633579624419452, 0.19021122030046342, -0.046485637341802263, -0.17641871952137458, 0.39833498112240728, 0.39523467764901271, -0.13364434540817743, -0.0011828142103154379, 0.18333362489004065, 0.072339887446854795, 0.096811654731024183, 0.2623225706024962, -0.038880711812973333, -0.095147356389043078, -0.32797689561579313, 0.38345998517535346, 0.40148749262155331, -0.39277937313364725, -0.40774631458499094, 0.36621151852507616, -0.40822442698887101, -0.23351079632264593, 0.33166073290449449, 0.36898082423156148, 0.37528853989229183, 0.24991809930911463, -0.3465346252896091, 0.016762855515073033, -0.078575213920283271, -0.31231181342519654, -0.40260968178095463, -0.1218748603257475, 0.31071908551606819, 0.1562931653949482, 0.10230613034434582, 0.127941328512514, -0.28492166302832456, 0.38966444184637822, 0.40769492078000907, 0.10528398564149932, -0.40206800725525388, -0.37636514267750426};
			valsToReplaceReal.insert(valsToReplaceReal.end(), toAdd.begin(), toAdd.begin()+d*d);
			valsToReplaceImag.insert(valsToReplaceImag.end(), toAdd.begin()+d*d, toAdd.end());
		}
	}

	// Combine the real and imag vecs
	std::vector<double> valsToReplace = {};
	valsToReplace.insert(valsToReplace.end(), valsToReplaceReal.begin(), valsToReplaceReal.end());
	valsToReplace.insert(valsToReplace.end(), valsToReplaceImag.begin(), valsToReplaceImag.end());

	// The indices to change
	for (int i=0; i<valsToReplaceReal.size(); i++) {
		indsToReplace.push_back(i);
	}
	for (int i=0; i<valsToReplaceImag.size(); i++) {
		indsToReplace.push_back(i+numVarsNonConj);
	}

	// For each of the equations
	for (int i=0; i<eqns.size(); i++) {

		// Substitute values if possible
		Polynomial<double> newEqn = eqns[i].substitute(indsToReplace, valsToReplace);

		// If the result is non-empty
		if (newEqn.coeffs.size() > 0) {

			// Check if this equation is unique
			bool unique = true;
			for (int j=0; j<newEqns.size(); j++) {
				if (newEqns[j] == newEqn) {
					unique = false;
					break;
				}
			}

			// If it is, add it
			if (unique) {
				newEqns.push_back(newEqn);
				//std::cout << i << ") 0 = " << newEqn << std::endl;
			}

		}

	}

	// Combine these to create a single polynomial
	std::cout << "Creating single polynomial..." << std::endl;
	Polynomial<double> poly(numVars);
	for (int i=0; i<newEqns.size(); i++) {
		poly += newEqns[i]*newEqns[i];
	}

	// Make sure all vars start from zero
	std::cout << "Simplifying variable indices..." << std::endl;
	std::vector<int> varList = poly.getVars();
	std::vector<int> mapTo;
	for (int i=0; i<varList.size(); i++) {
		mapTo.push_back(i);
	}
	poly = poly.changeVars(varList, mapTo); 
	std::cout << "Final poly has " << mapTo.size() << " variables" << std::endl;
	//std::cout << poly << std::endl;

	std::cout << std::endl;
	for (int i=0; i<newEqns.size(); i++) {
		newEqns[i] = newEqns[i].changeVars(varList, mapTo);
		std::cout << newEqns[i] << std::endl;
	}
	std::cout << std::endl;
	return 0;

	// Find a root of this
	std::cout << "Attempting to find a root..." << std::endl;
	poly.numVars += 1;
	std::vector<double> x = poly.findRoot(poly.numVars-1, 0.5, 1e-10, 10000000);
	//std::vector<double> x = poly.findRoot(2, 0.01, 1e-10, 10000000);
	std::cout << "Testing this x = " << poly.eval(x) << std::endl;
	for (int i=0; i<newEqns.size(); i++) {
		std::cout << "Testing eqn = " << newEqns[i].eval(x) << std::endl;
	}

	// Get the complex relaxation of this TODO
	//std::cout << "Attempting to find a relaxed root..." << std::endl;
	//Polynomial<double> relaxed = poly.getComplexRelaxation();
	//int origVars = poly.numVars;
	//relaxed.numVars += 1;
	//std::vector<double> x2 = relaxed.findRoot(relaxed.numVars-1);
	//std::vector<std::complex<double>> xComp(x2.size());
	//for (int i=0; i<origVars; i++) {
		//xComp[i] = std::complex<double>(x2[i], x2[i+origVars]);
		//std::cout << xComp[i] << std::endl;
	//}
	//std::cout << "Testing this x = " << relaxed.eval(x2) << std::endl;
	//std::cout << "Testing this x = " << Polynomial<std::complex<double>>(poly).eval(xComp) << std::endl;
	//for (int i=0; i<newEqns.size(); i++) {
		//std::cout << "Testing eqn = " << Polynomial<std::complex<double>>(newEqns[i]).eval(xComp) << std::endl;
	//}

	// Maybe a Groebner basis? TODO

	return 0;

}
	
