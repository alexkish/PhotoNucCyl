/////////////////////////////////////////////////////////
//													   //
//   ++++++++++++++++++++++++++++++++++++++++++++++++  //
//   + Deuterium Fusion Neutron Generator at UZH       //
//   + Alexander Kish                  UZH, 2008	   //
// 	 ++++++++++++++++++++++++++++++++++++++++++++++++  //
//													   //
/////////////////////////////////////////////////////////

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Orb.hh"
#include "G4Polyhedra.hh"
#include "G4Trd.hh"
#include "G4Cons.hh"

#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"

#include "G4VisAttributes.hh"

#include "G4OpBoundaryProcess.hh"
#include "G4SDManager.hh"

#include "G4Colour.hh"
#include "globals.hh"

#include "G4ios.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"

#include "vector"
#include "numeric"
#include "sstream"
#include "algorithm"
#include "cmath"
#include "cassert"

using std::vector;
using std::stringstream;
using std::max;

#include "NeutronLabConstruction.hh"
#include "ScreenSensitiveDetector.hh"


NeutronLabConstruction::NeutronLabConstruction()
{
}

NeutronLabConstruction::~NeutronLabConstruction()
{
}

G4VPhysicalVolume *NeutronLabConstruction::Construct()
{
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
G4cout 															<< G4endl;
G4cout << 	" ============================================= "	<< G4endl;
G4cout <<	"|   Photonuclear neutron source optimization  |"	<< G4endl;
G4cout <<	"|  ------------------------------------------ |"	<< G4endl;
G4cout <<	"|              Alexander Kish, UHM/CERN 2020  |"	<< G4endl;
G4cout <<	" ============================================= "	<< G4endl;
G4cout <<	"| "												<< G4endl;
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

	G4double	density,	// density
    			a,			// atomic weight (g/mole)
    			z;			// atomic number
	G4String	name,		// name
 				symbol;		// symbol
	G4int		ncomponents,
 				natoms;          
	G4double	temperature,
				pressure;    

	//----------	Define elements		--------------------------------------------------------
	G4Element *H  = new G4Element(name= "Hydrogen",		symbol= "H",  z= 1.,  a= 1.008   *CLHEP::g/CLHEP::mole);
	G4Element *B  = new G4Element(name= "Boron",		symbol= "B",  z= 5.,  a= 10.811  *CLHEP::g/CLHEP::mole);
	G4Element *C  = new G4Element(name= "Carbon",		symbol= "C",  z= 6.,  a= 12.011  *CLHEP::g/CLHEP::mole);
	G4Element *O  = new G4Element(name= "Oxygen",  		symbol= "O",  z= 8.,  a= 16.00   *CLHEP::g/CLHEP::mole);
	G4Element *F  = new G4Element(name= "Fluorine",		symbol= "F",  z= 9.,  a= 18.998  *CLHEP::g/CLHEP::mole);
	G4Element* Mg = new G4Element(name= "Magnesium", 	symbol= "Mg", z= 12., a= 24.305  *CLHEP::g/CLHEP::mole);
	G4Element *Al = new G4Element(name= "Aluminium", 	symbol= "Al", z= 13., a= 26.982  *CLHEP::g/CLHEP::mole);
	G4Element *Si = new G4Element(name= "Silicon",   	symbol= "Si", z= 14., a= 28.086  *CLHEP::g/CLHEP::mole);
	G4Element *P  = new G4Element(name= "Phosphorus",  	symbol= "P",  z= 15., a= 30.9738 *CLHEP::g/CLHEP::mole);
	G4Element *S  = new G4Element(name= "Sulphur",   	symbol= "S",  z= 16., a= 32.065  *CLHEP::g/CLHEP::mole);
	G4Element *K  = new G4Element(name= "Potassium", 	symbol= "K",  z= 19., a= 39.0983 *CLHEP::g/CLHEP::mole);
	G4Element* Ca = new G4Element(name= "Calcium",   	symbol= "Ca", z= 20., a= 40.078  *CLHEP::g/CLHEP::mole);
	G4Element *Ti = new G4Element(name= "Titanium",  	symbol= "Ti", z= 22., a= 47.867  *CLHEP::g/CLHEP::mole);
	G4Element *Cr = new G4Element(name= "Chromium",  	symbol= "Cr", z= 24., a= 51.996  *CLHEP::g/CLHEP::mole);
	G4Element *Mn = new G4Element(name= "Manganese", 	symbol= "Mn", z= 25., a= 54.938  *CLHEP::g/CLHEP::mole);
	G4Element *Fe = new G4Element(name= "Iron",			symbol= "Fe", z= 26., a= 55.845  *CLHEP::g/CLHEP::mole);
	G4Element *Co = new G4Element(name= "Cobalt",		symbol= "Co", z= 27., a= 58.9332 *CLHEP::g/CLHEP::mole);
	G4Element *Ni = new G4Element(name= "Nickel",    	symbol= "Ni", z= 28., a= 58.693  *CLHEP::g/CLHEP::mole);
	G4Element *Mo = new G4Element(name= "Molybdenum",	symbol= "Mo", z= 42., a= 95.94   *CLHEP::g/CLHEP::mole);
	G4Element *Pb = new G4Element(name= "Lead",			symbol= "Pb", z= 82., a= 207.2   *CLHEP::g/CLHEP::mole);
   	//G4Element* Ta  = new G4Element(name="Tallium",   	symbol="Ta", z=73, a=180.94*g/mole  );
   	//G4Element* W   = new G4Element(name="Tungsten",  	symbol="W",  z=74, a=182.30*g/mole  );
	G4Element* Li = G4NistManager::Instance()->FindOrBuildElement("Li");
	G4Element* Be 	= G4NistManager::Instance()->FindOrBuildElement("Be");
	G4Element* Ta 	= G4NistManager::Instance()->FindOrBuildElement("Ta");
	G4Element* W 	= G4NistManager::Instance()->FindOrBuildElement("W");

	
	//----------	 Define Materials 	--------------------------------------------------------

    // Air
	G4Material *Air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
	G4Material *N2 	= G4NistManager::Instance()->FindOrBuildMaterial("G4_NITROGEN");

    // Vacuum
	G4Material *Vacuum = new G4Material(name = "Vacuum", z = 1., a = 1. *CLHEP::g/CLHEP::mole, density = 1.e-20 *CLHEP::g/CLHEP::cm3, kStateGas, temperature = 0.1 * CLHEP::kelvin, pressure = 1.e-20 * CLHEP::bar);

    // Lead
    G4Material *Lead_mat = new G4Material(name = "Lead_mat", density = 11.34 *CLHEP::g/CLHEP::cm3, ncomponents =1);
    Lead_mat-> AddElement(Pb, 1); 

	// Concrete
	G4Material *Concrete = G4NistManager::Instance()->FindOrBuildMaterial("G4_CONCRETE");

	// Paraffin
	G4Material *Paraffin = G4NistManager::Instance()->FindOrBuildMaterial("G4_PARAFFIN");

    // Steel
	G4Material *Steel = new G4Material(name = "Steel", density = 7.7 *CLHEP::g/CLHEP::cm3, ncomponents = 3);
	Steel->AddElement(C,  0.04);
	Steel->AddElement(Fe, 0.88);
	Steel->AddElement(Co, 0.08);

	//------------------------------- stainless steel -------------------------------
	G4Material *Steel304 = new G4Material("Steel304", density = 8.00 *CLHEP::g/CLHEP::cm3, ncomponents =10);
	Steel304->AddElement(C,  0.0008);
	Steel304->AddElement(Si, 0.01);
	Steel304->AddElement(Mn, 0.02);
	Steel304->AddElement(P,  0.00045);
	Steel304->AddElement(S,  0.0003);
	Steel304->AddElement(Ni, 0.12);
	Steel304->AddElement(Cr, 0.17);
	Steel304->AddElement(Mo, 0.025);
	Steel304->AddElement(Ti, 0.004);
	Steel304->AddElement(Fe, 0.64945);

	//http://www.shieldwerx.com/poly-neutron.html
	G4Material *boron 	= G4NistManager::Instance()->FindOrBuildMaterial("G4_B"); // pure boron
	G4Material *PE 		= G4NistManager::Instance()->FindOrBuildMaterial("G4_POLYETHYLENE"); // pure polyethylene, e.g. SWX213
	G4Material *bismuth = G4NistManager::Instance()->FindOrBuildMaterial("G4_BISMUTH"); // pure bismuth
	// 5% borated polyethylene SWX201
	G4Material *BPE5 = new G4Material(name="BPE5", density=0.95*CLHEP::g/CLHEP::cm3, ncomponents=2); 
	BPE5->AddMaterial(boron, 	0.05);
	BPE5->AddMaterial(PE, 		0.95);
	// 5% borated high-density polyethylene SWX201HD
	G4Material *BPE5HD = new G4Material(name="BPE5HD", density=1.07*CLHEP::g/CLHEP::cm3, ncomponents=2); 
	BPE5HD->AddMaterial(boron, 	0.05);
	BPE5HD->AddMaterial(PE, 	0.95);
	// 30% borated polyethylene SWX210
	G4Material *BPE30 = new G4Material(name="BPE30", density=1.19*CLHEP::g/CLHEP::cm3, ncomponents=2); 
	BPE30->AddMaterial(boron, 	0.30);
	BPE30->AddMaterial(PE, 		0.70);
	// bismuth-based PE, poly-biz SWX217
/*	G4Material *BiPE = new G4Material(name="BiPE", density=2.92*CLHEP::g/CLHEP::cm3, ncomponents=2); 
	BiPE->AddMaterial(bismuth, 	0.215);
	BiPE->AddMaterial(PE, 		0.785);
*/    // Boric Acid
	G4Material *BoricAcid = new G4Material(name="BoricAcid", density=1.435*CLHEP::g/CLHEP::cm3, ncomponents =3);
	BoricAcid->AddElement(H, natoms=3);
	BoricAcid->AddElement(B, natoms=1);
	BoricAcid->AddElement(O, natoms=3);

	//dirt (LNGS rock for now)
	G4Material* rock = new G4Material("rock", 2.71*CLHEP::g/CLHEP::cm3, 8, kStateSolid);
	rock->AddElement(O,  0.5077);
	rock->AddElement(Ca, 0.2689);
	rock->AddElement(C,  0.1217);
	rock->AddElement(Mg, 0.0832);
	rock->AddElement(Si, 0.0105);
	rock->AddElement(Al, 0.0063);
	rock->AddElement(K,  0.0010);
	rock->AddElement(H,  0.0007);

	//Lithium Polyethylene
	//G4Material* polyethylene = man->FindOrBuildMaterial("G4_POLYETHYLENE");
	G4Material *LiPE = new G4Material("LiPE", 1.06*CLHEP::g/CLHEP::cm3, 2);
	LiPE->AddMaterial(PE, 92.46*CLHEP::perCent);
	LiPE->AddElement(Li, 7.54*CLHEP::perCent);

	// Li-6 thermal neutron absorber
	G4Isotope *Li6 = new G4Isotope("Lithium6", 3, 6, 6.015123*CLHEP::g/CLHEP::mole);
	G4Element *LiPure = new G4Element("LithiumPure", "LiPure", 1);
	LiPure->AddIsotope(Li6,1);
	G4Material *LiMat = new G4Material(name = "LiMat", density=0.534*CLHEP::g/CLHEP::cm3, ncomponents=1);
	LiMat->AddElement(LiPure,1);

	// S-32 neutron filter
	G4Isotope *S32 = new G4Isotope("Sulphur32", 16, 32, 31.972071*CLHEP::g/CLHEP::mole);
	G4Element *SPure = new G4Element("SulphurPure", "SiPure", 1);
	SPure->AddIsotope(S32,1);
	G4Material *SMat = new G4Material(name = "SMat", density=2.0*CLHEP::g/CLHEP::cm3, ncomponents=1);
	SMat->AddElement(SPure,1);

	// Pure silicon for neutron moderator material
 	G4Material *SiMat =  new G4Material(name="SiMat", density=2.33*CLHEP::g/CLHEP::cm3, ncomponents=1);
  	SiMat->AddElement(Si, 1);

	// Beryllium oxide for YBe source
	G4Material *BeO = new G4Material(name="BeO", density=3.01*CLHEP::g/CLHEP::cm3, ncomponents=2);
	BeO->AddElement(Be, 1);
	BeO->AddElement(O,  1);

  	// Pure beryllium
	G4Material *BeMat = new G4Material(name="BeMat", density=1.85*CLHEP::g/CLHEP::cm3, ncomponents=1);
	BeMat->AddElement(Be, 1.0);
  	// tungsten
	G4Material *WMat = new G4Material(name="WMat", density=19.25*CLHEP::g/CLHEP::cm3, ncomponents=1);
	WMat->AddElement(W, 1.0);

  	// tantalum
	G4Material *TaMat = new G4Material(name="TaMat", density=16.69*CLHEP::g/CLHEP::cm3, ncomponents=1);
	TaMat->AddElement(Ta, 1.0);



	////////////////////////////////////////////////////////////////
	//--------------	colours		--------------------------------
	G4Colour white  (1.0,	1.0,	1.0);
	G4Colour grey   (0.5,	0.5,	0.5);
	G4Colour lgrey  (.85,	.85,	.85);
	G4Colour red    (1.0,	0.0,	0.0);
	G4Colour lred   (0.75,	0.0,	0.0);
	G4Colour xlred  (0.5,	0.0,	0.0);
	G4Colour cyan   (0.0,	1.0,	1.0);
	G4Colour blue   (0.0,	0.0,	1.0);
	G4Colour lblue  (.5,	0.5,	1.,		1.);
	G4Colour xlblue (.5,	0.5,	1.,		0.2);
	G4Colour magenta(1.0,	0.0,	1.0);
	G4Colour yellow (1.0,	1.0,	0.0);
	G4Colour green  (0.,	.1,		0.);
	G4Colour lgreen (0.0,	.75,	0.0);
	G4Colour xlgreen(0.0,	0.5,	0.0);
	G4Colour brown  (0.7,	0.4,	0.1);
	G4Colour orange (1.0,	0.5,	0.0);
	G4Colour xlorange (1.0,	0.5,	0.0, 	0.2);

	//--------------	rotations	--------------------------------
	G4RotationMatrix *RotationXPlus225	= new G4RotationMatrix();
	RotationXPlus225->rotateX(22.5*CLHEP::deg);

	G4RotationMatrix *RotationXMinus225	= new G4RotationMatrix();
	RotationXMinus225->rotateX(-22.5*CLHEP::deg);

	G4RotationMatrix *RotationXPlus45	= new G4RotationMatrix();
	RotationXPlus45->rotateX(45.*CLHEP::deg);

	G4RotationMatrix *RotationXMinus45	= new G4RotationMatrix();
	RotationXMinus45->rotateX(-45.*CLHEP::deg);

	G4RotationMatrix *RotationXPlus90	= new G4RotationMatrix();
	RotationXPlus90->rotateX(90.*CLHEP::deg);
	RotationXPlus90->rotateY(0.*CLHEP::deg);
	RotationXPlus90->rotateZ(0.*CLHEP::deg);

	G4RotationMatrix *RotationYPlus90	= new G4RotationMatrix();
	RotationYPlus90->rotateX(0.*CLHEP::deg);
	RotationYPlus90->rotateY(90.*CLHEP::deg);
	RotationYPlus90->rotateZ(0.*CLHEP::deg);

	G4RotationMatrix *RotationXMinus90	= new G4RotationMatrix();
	RotationXMinus90->rotateX(-90.*CLHEP::deg);

	G4RotationMatrix *RotationX180		= new G4RotationMatrix();
	RotationX180->rotateX(180.*CLHEP::deg);
	
//	G4RotationMatrix *ZeroRot			= new G4RotationMatrix();
//	ZeroRot->rotateX(0. *deg);
//	ZeroRot->rotateY(0. *deg);
//	ZeroRot->rotateZ(0. *deg);
	
	G4RotationMatrix ZeroRot;
	ZeroRot.rotateX(0. *CLHEP::deg);
	ZeroRot.rotateY(0. *CLHEP::deg);
	ZeroRot.rotateZ(0. *CLHEP::deg);


	G4double opendeg	= 0.0 *CLHEP::deg;
	G4double closedeg	= 360.0 *CLHEP::deg;

	////////////////////////////////////////
	// lab hall
	G4double Lab_HalfX = 1.5 *CLHEP::m;
	G4double Lab_HalfY = 1.5 *CLHEP::m;
	G4double Lab_HalfZ = 1.5 *CLHEP::m;

	G4Box *Laboratory = new G4Box("Laboratory", Lab_HalfX, Lab_HalfY, Lab_HalfZ);

	G4double pipe_thick 		= 1.5 *CLHEP::mm;
	G4double pipe_OD 			= 3.0 *CLHEP::cm;
	G4double pipe_ID 			= 0;
	G4double pipe_height 		= 10.0 *CLHEP::cm;

	// shield for internal source
	G4double shield_thick 		= 10.0 * CLHEP::cm;
	G4double shield_ID 			= 0;
	G4double shield_OD 			= pipe_OD + 2*shield_thick;
	G4double shield_height 		= pipe_height;

	G4double interior_OD 		= pipe_OD - 2*pipe_thick;
	G4double interior_ID 		= 0;
	G4double interior_height 	= pipe_height;

	G4double core_ID 			= 0;
	G4double core_OD 			= 1.0 *CLHEP::mm;

	// sphere for interior use
	G4double shell_thick 		= 1.0 * CLHEP::cm;
	G4double shell_ID 			= 0;
	G4double shell_OD 			= core_OD + 2*shell_thick;

	// cylinder for exterior use
	G4double shellTube_thick 	= 2.0 * CLHEP::cm;
	G4double shellTube_ID 		= 0;
	G4double shellTube_OD 		= pipe_OD + 2*shellTube_thick;
	G4double shellTube_height	= pipe_height;

/*	// shield for external source
	G4double shield_thick 		= 10.0 * CLHEP::cm;
	G4double shield_ID 			= 0;
	G4double shield_OD 			= shellTube_OD + 2*shield_thick;
	G4double shield_height 		= pipe_height;
*/
	G4double screen_thick 		= 0.001 *CLHEP::mm;
	G4double screen_ID 			= 0;
	G4double screen_OD 			= shield_OD + 2*screen_thick;
	G4double screen_height 		= pipe_height;

	G4Tubs *screen		= new G4Tubs("screen", screen_ID/2, screen_OD/2, screen_height/2, opendeg, closedeg);	
	G4Tubs *shield		= new G4Tubs("shield", shield_ID/2, shield_OD/2, shield_height/2, opendeg, closedeg);	
	G4Tubs *pipe		= new G4Tubs("pipe", pipe_ID/2, pipe_OD/2, pipe_height/2, opendeg, closedeg);	
	G4Tubs *interior	= new G4Tubs("interior", interior_ID/2, interior_OD/2, interior_height/2, opendeg, closedeg);	
	G4Orb  *core 		= new G4Orb("core", core_OD/2);
	G4Orb  *shell 		= new G4Orb("shell", shell_OD/2);
	G4Tubs *shellTube	= new G4Tubs("shellTube", shellTube_ID/2, shellTube_OD/2, shellTube_height/2, opendeg, closedeg);	

	// screen to detect particles

	///////////////////////////////////////////////////////////////
	// Logical Volumes (declared in 'NeutronLabDetectorGeometry.hh')
	G4LogicalVolume *Laboratory_log 	= new G4LogicalVolume(Laboratory, 	Air, 	"Laboratory_log");
	G4LogicalVolume *screen_log 		= new G4LogicalVolume(screen, 		Air, 	"screen_log");
	G4LogicalVolume *shield_log 		= new G4LogicalVolume(shield, 		WMat, 	"shield_log");
	G4LogicalVolume *pipe_log 			= new G4LogicalVolume(pipe, 		Steel, 	"pipe_log");
	G4LogicalVolume *interior_log 		= new G4LogicalVolume(interior, 	Air, 	"interior_log");
	G4LogicalVolume *shell_log 			= new G4LogicalVolume(shell, 		BeO, 	"shell_log");
	G4LogicalVolume *shellTube_log 		= new G4LogicalVolume(shellTube, 	BeO, 	"shellTube_log");
	G4LogicalVolume *core_log 			= new G4LogicalVolume(core, 		Steel, 	"core_log");

	// Physical Volumes
	// internal source
	G4PVPlacement *Laboratory_phys 		= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), Laboratory_log, "Laboratory_phys", 	0, false, 0);
	G4PVPlacement *screen_phys 			= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), screen_log, 		"screen_phys", 		Laboratory_log, false, 0);
	G4PVPlacement *shield_phys 			= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), shield_log, 		"shield_phys", 		screen_log, false, 0);
	G4PVPlacement *pipe_phys 			= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), pipe_log, 		"pipe_phys", 		shield_log, false, 0);
	G4PVPlacement *interior_phys 		= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), interior_log, 	"interior_phys", 	pipe_log, false, 0);
	G4PVPlacement *shell_phys 			= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), shell_log, 		"shell_phys", 		interior_log, false, 0);
	G4PVPlacement *core_phys 			= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), core_log, 		"core_phys", 		shell_log, false, 0);

/*	// external source
	G4PVPlacement *Laboratory_phys 		= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), Laboratory_log, "Laboratory_phys", 	0, false, 0);
	G4PVPlacement *screen_phys 			= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), screen_log, 		"screen_phys", 		Laboratory_log, false, 0);
	G4PVPlacement *shield_phys 			= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), shield_log, 		"shield_phys", 		screen_log, false, 0);
	G4PVPlacement *shellTube_phys 		= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), shellTube_log, 	"shellTube_phys", 	shield_log, false, 0);
	G4PVPlacement *pipe_phys 			= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), pipe_log, 		"pipe_phys", 		shellTube_log, false, 0);
	G4PVPlacement *interior_phys 		= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), interior_log, 	"interior_phys", 	pipe_log, false, 0);
	G4PVPlacement *core_phys 			= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), core_log, 		"core_phys", 		interior_log, false, 0);
*/

	///////////////////////////////////////////////////////////////////////
	// Visual Attributes
	Laboratory_log				->SetVisAttributes(G4VisAttributes::Invisible);		

	G4VisAttributes *screen_vis = new G4VisAttributes(xlgreen);
	screen_vis					->SetVisibility(true);
	screen_vis					->SetForceSolid(true);
	screen_log					->SetVisAttributes(screen_vis);

	G4VisAttributes *shield_vis = new G4VisAttributes(grey);
	shield_vis					->SetVisibility(true);
	shield_vis					->SetForceSolid(true);
	shield_log					->SetVisAttributes(shield_vis);

	G4VisAttributes *pipe_vis = new G4VisAttributes(blue);
	pipe_vis					->SetVisibility(true);
	pipe_vis					->SetForceSolid(true);
	pipe_log					->SetVisAttributes(pipe_vis);

	G4VisAttributes *shell_vis = new G4VisAttributes(orange);
	shell_vis					->SetVisibility(true);
	shell_vis					->SetForceSolid(true);
	shell_log					->SetVisAttributes(shell_vis);
	shellTube_log				->SetVisAttributes(shell_vis);

	G4VisAttributes *core_vis = new G4VisAttributes(red);
	core_vis					->SetVisibility(true);
	core_vis					->SetForceSolid(true);
	core_log					->SetVisAttributes(core_vis);

	// get mass of the shell
	G4double mass_shell = 0.;
	mass_shell += shell_log ->GetMass(false, false)/CLHEP::kg;
	G4cout 	<< "Mass of the shell :" << mass_shell << " kg" << G4endl;

	// get mass of the shell tube
	G4double mass_shellTube = 0.;
	mass_shellTube += shellTube_log ->GetMass(false, false)/CLHEP::kg;
	G4cout 	<< "Shell thick:            " << shellTube_thick/CLHEP::cm << " cm" << G4endl;
	G4cout 	<< "Mass of the shell tube: " << mass_shellTube << " kg" << G4endl;

	// get mass of the shield
	G4double mass_shield = 0.;
	mass_shield += shield_log ->GetMass(false, false)/CLHEP::kg;
	G4cout 	<< "Shield thick:       " << shield_thick/CLHEP::cm << " cm" << G4endl;
	G4cout 	<< "Shield OD:          " << shield_OD/CLHEP::cm << " cm" << G4endl;
	G4cout 	<< "Mass of the shield: " << mass_shield << " kg" << G4endl;



	// SENSITIVE DETECTOR
	G4SDManager *pSDManager = G4SDManager::GetSDMpointer();

	// make entire laboratory sensitive
	ScreenSensitiveDetector *pScreenSD = new ScreenSensitiveDetector("ScreenSD");
	pSDManager->AddNewDetector(pScreenSD);
	screen_log	->SetSensitiveDetector(pScreenSD);


//--------------	END PROGRAM	------------------------
return Laboratory_phys;
//------------------------------------------------------

}


//NeutronLabDetectorConstruction::~NeutronLabDetectorConstruction()
//{
//delete m_pDetectorMessenger;
//}








