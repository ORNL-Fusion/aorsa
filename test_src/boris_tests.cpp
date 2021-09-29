#include <iostream>
#include <thrust/execution_policy.h>
#include "test_utils.hpp"
#include "config_interface.h"
#include "test_data_filepath.hpp"
#include "utils.h"
#include "flags.hpp"
#include "Particles.h"
#include "boris.h"
#include "Surfaces.h"
#include "geometryCheck.h"
#include "boundaryInit.h"

template <typename T=double>
bool compareVectors(std::vector<T> a, std::vector<T> b, T epsilon, T margin)
{
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); i++) 
  {
    
    bool margin_check = (a[i] != Approx(b[i]).margin(margin));
    bool epsilon_check = (a[i] != Approx(b[i]).epsilon(epsilon));
    
    if (margin_check && epsilon_check)
    {
      
      std::cout << "margin epsilon " <<
        margin_check << " " << epsilon_check << std::endl; 
      std::cout << "Element " << i << 
        " " << a[i] << " Should == " << b[i] << std::endl;
      
      return false;
    }
  }
  
  return true;
}

TEST_CASE( "boris - not fully implemented" )
{
  SECTION( "e cross b" )
  {
    /* timesteps */
    int nT = 1e6;
    //int nT = 1;

    /* create particles */
    libconfig::Config cfg_geom;

    cfg_geom.setAutoConvert(true);

    importLibConfig(cfg_geom, BORIS_TEST_FILE);

    auto gitr_flags = new Flags( cfg_geom );

    libconfig::Setting &impurity = cfg_geom.lookup( "impurityParticleSource" );

    /* 2nd argument is deprecated - random number stream related */
    auto particleArray = new Particles( 1, 1, cfg_geom, gitr_flags );

    /* Captain! Important equation to convert eV energy to vector velocity components */
    gitr_precision E = 14;
    gitr_precision amu = 27;
    gitr_precision vtotal = std::sqrt(2.0 * E * 1.602e-19 / amu / 1.66e-27);
    std::cout << "vtotal: " << vtotal << std::endl;

    /* set particle properties: */
    gitr_precision dt = 1.0e-6;
    particleArray->setParticleV( 0, 0, 0, 0, vtotal, 0, 0, 13, amu, 2.0, dt );
    /*
  void setParticleV(int indx, gitr_precision x, gitr_precision y, gitr_precision z,
                    gitr_precision Vx, gitr_precision Vy, gitr_precision Vz,
                    gitr_precision Z, gitr_precision amu, gitr_precision charge,
                    gitr_precision dt)
    */

    thrust::counting_iterator<std::size_t> particle_iterator_start(0);
    thrust::counting_iterator<std::size_t> particle_iterator_end(1);

    /* Captain... Just rip everything from cross field diffusion tests and convert it
       into a boris test. Compare and contrast differences between what's in gitr.cpp
       and what's in cross_field_diffusion_tests.cpp */
    int nLines = 2;
    sim::Array<Boundary> boundaries( nLines + 1, Boundary() );

    int nSurfaces = importGeometry( cfg_geom, boundaries );

    REQUIRE( nSurfaces == 2 );

    /* start */
    int nHashes = 1;
    sim::Array<int> nR_closeGeom(nHashes, 0);
    sim::Array<int> nY_closeGeom(nHashes, 0);
    sim::Array<int> nZ_closeGeom(nHashes, 0);
    sim::Array<int> nHashPoints(nHashes, 0);
    sim::Array<int> n_closeGeomElements(nHashes, 0);
    int nEdist = 1;
    gitr_precision E0dist = 0.0;
    gitr_precision Edist = 0.0;
    int nAdist = 1;
    gitr_precision A0dist = 0.0;
    gitr_precision Adist = 0.0;

    auto surfaces = new Surfaces(nSurfaces, nEdist, nAdist);
    surfaces->setSurface(nEdist, E0dist, Edist, nAdist, A0dist, Adist);
    sim::Array<gitr_precision> closeGeomGridr(1),
      closeGeomGridy(1), closeGeomGridz(1);
    sim::Array<int> closeGeom(1, 0);
    /* end */
    /* 
       USE_ADAPTIVE_DT = 0 run with both options, should get same answer. Start with 0 
       turn off all hashing stuff
    */

    gitr_precision perpDiffusionCoeff = 0.0;
    cfg_geom.lookupValue("backgroundPlasmaProfiles.Diffusion.Dperp",
                          perpDiffusionCoeff);

    int nR_Bfield = 1, nZ_Bfield = 1, n_Bfield = 1;

    /* required option: USE_PRESHEATH_EFIELD=1 and GITR_BFIELD_INTERP=1 */
    /* create a unified setup script */
    sim::Array<gitr_precision> br(n_Bfield), by(n_Bfield), bz(n_Bfield);

    /* uniform bfield */
    br[ 0 ] = 0;
    /* large bfield in teslas gives smaller gyromotion radius */
    by[ 0 ] = 5;
    bz[ 0 ] = 0;

    /* for the uniform efield, set efield to 1000 in z just make the cross product geometry */
    /* presheath efield is in the bulk plasma and sheath efield is at the surface of the wall */

    sim::Array<gitr_precision> bfieldGridr(nR_Bfield), bfieldGridz(nZ_Bfield);

    int nR_PreSheathEfield = 1;
    int nY_PreSheathEfield = 1;
    int nZ_PreSheathEfield = 1;

    int nPSEs = nR_PreSheathEfield * nY_PreSheathEfield * nZ_PreSheathEfield;

    sim::Array<gitr_precision> preSheathEGridr(nR_PreSheathEfield),
      preSheathEGridy(nY_PreSheathEfield), preSheathEGridz(nZ_PreSheathEfield);

    sim::Array<gitr_precision> PSEr(nPSEs), PSEz(nPSEs), PSEt(nPSEs);
    PSEr[ 0 ] = 0;
    PSEz[ 0 ] = -1000;
    /* y and t */
    PSEt[ 0 ] = 0;

    int nR_closeGeom_sheath = 1;
    int nY_closeGeom_sheath = 1;
    int nZ_closeGeom_sheath = 1;
    int n_closeGeomElements_sheath = 1;
    int nGeomHash_sheath = 1;
    sim::Array<gitr_precision> closeGeomGridr_sheath(nR_closeGeom_sheath),
      closeGeomGridy_sheath(nY_closeGeom_sheath),
      closeGeomGridz_sheath(nZ_closeGeom_sheath);
    sim::Array<int> closeGeom_sheath(nGeomHash_sheath);

    /* create boris operator */
    move_boris boris( particleArray, dt, boundaries.data(), nLines, nR_Bfield, nZ_Bfield,
                      bfieldGridr.data(), bfieldGridz.data(), br.data(), bz.data(), by.data(),
                      nR_PreSheathEfield, 
                      nY_PreSheathEfield,
                      nZ_PreSheathEfield,
                      &preSheathEGridr.front(), &preSheathEGridy.front(),
                      &preSheathEGridz.front(), &PSEr.front(), &PSEz.front(), &PSEt.front(),
                      nR_closeGeom_sheath, nY_closeGeom_sheath, nZ_closeGeom_sheath,
                      n_closeGeomElements_sheath, closeGeomGridr_sheath.data(),
                      &closeGeomGridy_sheath.front(), &closeGeomGridz_sheath.front(),
                      &closeGeom_sheath.front(), gitr_flags );

    geometry_check geometry_check0(
        particleArray, nLines, &boundaries[0], surfaces, dt, nHashes,
        nR_closeGeom.data(), nY_closeGeom.data(), nZ_closeGeom.data(),
        n_closeGeomElements.data(), &closeGeomGridr.front(),
        &closeGeomGridy.front(), &closeGeomGridz.front(), &closeGeom.front(),
        nEdist, E0dist, Edist, nAdist, A0dist, Adist);

    /* get particle xyz before */
    /* time loop */
    std::cout << "Captain! num particles: " << particleArray->nParticles << std::endl;
    std::cout << "Captain! Before: " << particleArray->x[0] << " " << particleArray->z[0]
              << " " << particleArray->y[0]
              << std::endl;

    for (int tt = 0; tt < nT; tt++)
    {

      thrust::for_each( thrust::device,
                        particle_iterator_start,
                        particle_iterator_end,
                        boris );

      thrust::for_each(thrust::device,
                       particle_iterator_start, 
                       particle_iterator_end, 
                       geometry_check0 );
    }

    std::cout << "Captain! After: " << particleArray->x[0] << " " << particleArray->z[0]
              << " " << particleArray->y[0]
              << std::endl;
    /* get particle xyz after */

    /* what should it be analytically? */

    /* charge coulombs, e in V/m, b in Teslas */
    /* q * e / b^2 */
  }
  SECTION( "getE tests" )
  {
    libconfig::Config cfg_geom;

    cfg_geom.setAutoConvert(true);

    importLibConfig(cfg_geom, "../test_data/getE.cfg");
    int nLines = 1;
    sim::Array<Boundary> boundaries( nLines + 1, Boundary() );

    int nSurfaces = importGeometry( cfg_geom, boundaries );
    int nR_Dens = 1;
    int nZ_Dens = 1;
    sim::Array<gitr_precision> DensGridr(1, 0.0);
    sim::Array<gitr_precision> DensGridz(1, 0.0);
    sim::Array<gitr_precision> ni(1, 1.0e19);
    sim::Array<gitr_precision> ne(1, 1.0e19);
    
    // Temperature = 20 eV
    int nR_Temp = 1;
    int nZ_Temp = 1;
    sim::Array<gitr_precision> TempGridr(1, 0.0);
    sim::Array<gitr_precision> TempGridz(1, 0.0);
    sim::Array<gitr_precision> ti(1,20.0);
    sim::Array<gitr_precision> te(1,20.0);
    
    int nR_Bfield = 1, nZ_Bfield = 1, n_Bfield = 1;

    /* required option: USE_PRESHEATH_EFIELD=1 and GITR_BFIELD_INTERP=1 */
    /* create a unified setup script */
    sim::Array<gitr_precision> br(n_Bfield), by(n_Bfield), bz(n_Bfield);

    /* uniform bfield */
    br[ 0 ] = std::cos(M_PI*5.0/180);
    /* large bfield in teslas gives smaller gyromotion radius */
    by[ 0 ] = 0;
    bz[ 0 ] = -std::sin(M_PI*5.0/180);;

    /* for the uniform efield, set efield to 1000 in z just make the cross product geometry */
    /* presheath efield is in the bulk plasma and sheath efield is at the surface of the wall */

    sim::Array<gitr_precision> bfieldGridr(nR_Bfield), bfieldGridz(nZ_Bfield);
    gitr_precision background_Z = 1;
    gitr_precision background_amu = 2;
    gitr_precision biasPotential = 0;

  
    std::for_each(boundaries.begin(), boundaries.end() - 1,
                boundary_init(background_Z, background_amu, nR_Dens, nZ_Dens,
                              DensGridr.data(), DensGridz.data(), ni.data(),
                              ne.data(), nR_Bfield, nZ_Bfield,
                              bfieldGridr.data(), bfieldGridz.data(), br.data(),
                              bz.data(), by.data(), nR_Temp, nZ_Temp,
                              TempGridr.data(), TempGridz.data(), ti.data(),
                              te.data(), biasPotential));
    
    int nHashes = 1;
    int nR_closeGeom_sheath = 1;
    int nY_closeGeom_sheath = 1;
    int nZ_closeGeom_sheath = 1;
    int nHashPoints_sheath = 1;
    int n_closeGeomElements_sheath = 1;
    sim::Array<gitr_precision> closeGeomGridr_sheath(1),
      closeGeomGridy_sheath(1), closeGeomGridz_sheath(1);
    sim::Array<int> closeGeom_sheath(1, 0);
    
    int closestBoundaryIndex = 0;
    int surfIndex = 0;
    gitr_precision minDistance = 0.0;
    gitr_precision thisE[3] = {0.0};
    sim::Array<gitr_precision> px(1, 0);
    sim::Array<gitr_precision> py(1, 0);
    sim::Array<gitr_precision> pz(1, 0.001);
    gitr_precision dz = 0.005/10000.0;
     int nZ = 10000;
     std::vector<gitr_precision> gitrE(nZ,0.0);
    for(int j=0;j<nZ;j++)
    {
      pz[0] = j*dz;
      minDistance =
          getE(px[0], py[0], pz[0], thisE, boundaries.data(), nLines,
               nR_closeGeom_sheath, nY_closeGeom_sheath, nZ_closeGeom_sheath,
               n_closeGeomElements_sheath, &closeGeomGridr_sheath.front(),
               &closeGeomGridy_sheath.front(), &closeGeomGridz_sheath.front(),
               &closeGeom_sheath.front(), closestBoundaryIndex);
      gitrE[j] = thisE[2];
    }
      std::cout << "minDist " << minDistance << std::endl; 
      std::cout << "Efield " << thisE[0] << " " << thisE[1] << " " << thisE[2] << std::endl; 

std::ifstream in("../test_data/Efield.txt");

std::string str;
std::vector<gitr_precision> gold;
// Read the next line from File untill it reaches the end.
while (std::getline(in, str))
{
    // Line contains string of length > 0 then save it in vector
    if(str.size() > 0){
      std::size_t sz = str.length();
      gitr_precision val = std::atof(str.c_str());
        gold.push_back(-val);
        }
}
for(int i=0;i<gold.size();i++){
std::cout << gold[i] << std::endl;
}
gold[0] = 0.0;
    // Compare vectors to ensure reproducibility
    gitr_precision margin = 0.1;
    gitr_precision epsilon = 0.001;
    REQUIRE(compareVectors<gitr_precision>(gitrE,gold,epsilon,margin));
  }
}
