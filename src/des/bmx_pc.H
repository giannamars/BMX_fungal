//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#ifndef _BMX_PC_H_
#define _BMX_PC_H_

#include <AMReX_AmrParGDB.H>
#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>

#include <bmx_bc_list.H>
#ifdef NEW_CHEM
#include <bmx_chem.H>
#endif

/**
 * Note: if the enumeration starts at 0, then count will represent the total
 * number of values defined prior to the count variable.
 */

#ifndef NEW_CHEM
struct realData
{
  enum { // Real particle data stored in p.m_data
    velx,
    vely,
    velz,
    radius,
    volume,
    fluid_A,
    fluid_B,
    consume_A,
    consume_B,
    count
  };
};
#endif

#ifndef NEW_CHEM
struct intData
{
  enum { // Integer particle data stored in p.m_data
    state,
    phase,
    count
  };
};
#endif

///
/// We add a particle to the neighbor list if its center is within 3 particle radii
///
struct pairDebug
{
  enum {TIPFUSION, INTERIORFUSION, CLEANUP, IGNORE};
};

struct BMXCheckPair
{
  amrex::Real neighborhood;
  int locale;
  bool single_pass;
  /* neighborhood: square of cutoff distance
   * single_pass: only include pairs that fall within cutoff and for which i<j
   * locale: debugging flag to signal debugging code inside check pair routine
   */
  BMXCheckPair (amrex::Real neighborhood_) : neighborhood(neighborhood_), single_pass(true), locale(pairDebug::IGNORE) {}
  BMXCheckPair (amrex::Real neighborhood_, bool single_pass_) : neighborhood(neighborhood_),
      single_pass(single_pass_), locale(pairDebug::IGNORE) {}
  BMXCheckPair (amrex::Real neighborhood_, bool single_pass_, int locale_) : neighborhood(neighborhood_),
      single_pass(single_pass_), locale(locale_) {}

  template <class P>
  AMREX_GPU_DEVICE AMREX_INLINE
  bool operator() (const P* particles, const int i, const int j) const
  {
    const P& p1 = particles[i];
    const P& p2 = particles[j];
    int idx = p1.idata(intIdx::id);
    int icpu = p1.idata(intIdx::cpu);
    int jdx = p2.idata(intIdx::id);
    int jcpu = p2.idata(intIdx::cpu);

    if (single_pass) {
      if (i < j) {
        if ((p1.pos() - p2.pos()).radSquared() <= neighborhood)
          return true;  
      }
    } else { 
#if 0
      if (icpu == 6 && jcpu == 6) {
        if (idx == 318 || idx == 247 ||
            jdx == 247 || jdx == 318) {
          if (locale == pairDebug::TIPFUSION) {
          std::printf("Checking pair Tip Fusion idx: %d icpu: %d jdx: %d jcpu: %d r: %e nghb: %e\n",
              idx,icpu,jdx,jcpu,(p1.pos() - p2.pos()).radSquared(),neighborhood);
          } else if (locale == pairDebug::INTERIORFUSION) {
          std::printf("Checking pair Interior Fusion idx: %d icpu: %d jdx: %d jcpu: %d r: %e nghb: %e\n",
              idx,icpu,jdx,jcpu,(p1.pos() - p2.pos()).radSquared(),neighborhood);
          } else if (locale == pairDebug::CLEANUP) {
#endif
#if 0
          if (locale == pairDebug::CLEANUP) {
          std::printf("Checking pair Cleanup idx: %d icpu: %d inew: %d"
              " isplit: %d ifuse: %d jdx: %d jcpu: %d jnew: %d jsplit: %d"
              " jfuse: %d r: %e nghb: %e\n",
              idx,icpu,p1.idata(intIdx::new_flag),
              p1.idata(intIdx::split_flag),p1.idata(intIdx::fuse_flag),
              jdx,jcpu,p2.idata(intIdx::new_flag),
              p2.idata(intIdx::split_flag),p1.idata(intIdx::fuse_flag),
              (p1.pos() - p2.pos()).radSquared(),neighborhood);
          }
#endif
#if 0
          }
        }
      }
#endif
      if ((p1.pos() - p2.pos()).radSquared() <= neighborhood) {
        return true;  
      }
    }

    return false;
  }
};

struct BMXCheckPolyPair
{
  int typeIdx;
  int num_bin_types;
  amrex::Real *neighborhood;
  int locale;
  bool single_pass;
  /* neighborhood: square of cutoff distance
   * single_pass: only include pairs that fall within cutoff and for which i<j
   * locale: debugging flag to signal debugging code inside check pair routine
   */
  BMXCheckPolyPair (int typeIdx_, int num_bin_types_, amrex::Real *neighborhood_)
    : typeIdx(typeIdx_), num_bin_types(num_bin_types_),
      neighborhood(neighborhood_), single_pass(true),
      locale(pairDebug::IGNORE) {}
  BMXCheckPolyPair (int typeIdx_, int num_bin_types_, amrex::Real *neighborhood_,
      bool single_pass_) : typeIdx(typeIdx_), num_bin_types(num_bin_types_),
      neighborhood(neighborhood_), single_pass(single_pass_),
      locale(pairDebug::IGNORE) {}
  BMXCheckPolyPair (int typeIdx_, int num_bin_types_, amrex::Real *neighborhood_,
      bool single_pass_, int locale_) : typeIdx(typeIdx_),
      num_bin_types(num_bin_types_), neighborhood(neighborhood_),
      single_pass(single_pass_), locale(locale_) {}

  template <class P>
  AMREX_GPU_DEVICE AMREX_INLINE
  bool operator() (const P* particles, const int i, const int j, int type,
      bool ghost_i, bool ghost_j) const
  {
    const P& p1 = particles[i];
    const P& p2 = particles[j];
    int type_i = p1.idata(typeIdx);
    int type_j = p2.idata(typeIdx);
    if (!(ghost_i && ghost_j)) {    // One particle must be locally held

      if (single_pass) {
        if (i < j || type_i < type_j) { // Smaller index or type finds neighbor
          // Get index for neighborhood from analytical partial sums
          int imin = amrex::min(type_i,type_j);
          int imax = amrex::max(type_i,type_j);
          int k = imin-1;
          // \sum_{i=N-k}^{N}(i)
          int indx = num_bin_types*(k+1)-(k*(k+1))/2;
          // \sum_{i=1}^{l}(i)
          indx += imax-imin;
          if ((p1.pos() - p2.pos()).radSquared() <= neighborhood[indx])
            return true;
        }
      } else { 
        // Get index for neighborhood from analytical partial sums
        int imin = amrex::min(type_i,type_j);
        int imax = amrex::max(type_i,type_j);
        int k = imin-1;
        // \sum_{i=N-k}^{N}(i)
        int indx = num_bin_types*(k+1)-(k*(k+1))/2;
        // \sum_{i=1}^{l}(i)
        indx += imax-imin;
        if ((p1.pos() - p2.pos()).radSquared() <= neighborhood[indx])
          return true;  
      }
    }

    return false;
  }
};

// Cannot inherit from this class?
#ifdef NEW_CHEM
class BMXParticleContainer final
    : public amrex::NeighborParticleContainer<MAX_CHEM_REAL_VAR,MAX_CHEM_INT_VAR>
#else
class BMXParticleContainer final
    : public amrex::NeighborParticleContainer<realData::count,intData::count>
#endif
{
public:
#ifdef NEW_CHEM
    using BMXParIter = amrex::ParIter<MAX_CHEM_REAL_VAR,MAX_CHEM_INT_VAR,0,0>;
    using BMXParConstIter = amrex::ParConstIter<MAX_CHEM_REAL_VAR,MAX_CHEM_INT_VAR,0,0>;
#else
    using BMXParIter = amrex::ParIter<realData::count,intData::count,0,0>;
    using BMXParConstIter = amrex::ParConstIter<realData::count,intData::count,0,0>;
#endif

    // Get number of particles
    int NumberOfParticles (BMXParIter& pti)
    { return pti.numParticles(); }

    BMXParticleContainer (amrex::AmrCore* amr_core);

    virtual ~BMXParticleContainer () {};

    void InitData ();

    void InitParticlesAscii (const std::string& file);

    void WriteToAscii (const std::string& file,int nstep,amrex::Real time);
    void CountParticleTypes(int &nparticles, int &nsegments);
    void WriteSegmentsToAscii (const std::string& file,int nstep,
        amrex::Real time, int nsegments);
    void WriteParticlesToAscii (const std::string& file,int nstep,
        amrex::Real time, int nparticles);

    void printParticles ();

    Real computeParticleVolume () const;

    Real computeParticleContent (int comp) const;

    void AllocData ();

    void PrintParticleCounts ();

    void EvolveParticles (amrex::Real dt,
                          const amrex::Vector<amrex::MultiFab *> cost,
                          std::string & knapsack_weight_type,
                          int & nsubsteps);
    void ParticleExchange (Real dt,
                           const amrex::Vector<amrex::MultiFab*> cost,
                           std::string& knapsack_weight_type,
                           int& nsubsteps);

    bool EvaluateTipFusion(const amrex::Vector<amrex::MultiFab*> cost,
                                   std::string& knapsack_weight_type);

    void EvaluateInteriorFusion(const amrex::Vector<amrex::MultiFab*> cost,
                                   std::string& knapsack_weight_type);

    void CleanupFusion(const amrex::Vector<amrex::MultiFab*> cost,
                                   std::string& knapsack_weight_type);

    void InitBonds(const amrex::Vector<amrex::MultiFab*> cost,
                                   std::string& knapsack_weight_type);

    void PrintConnectivity(const amrex::Vector<amrex::MultiFab*> cost,
                                   std::string& knapsack_weight_type);

    void CalculateFungalCM(const amrex::Vector<amrex::MultiFab*> cost,
                                  std::string& knapsack_weight_type,
                                  amrex::RealVect &cm);

    void CalculateFungalRG(const amrex::Vector<amrex::MultiFab*> cost,
                                  std::string& knapsack_weight_type,
                                  Real &Rg, Real &Masst);
    void CalculateFungalDensityProfile(const amrex::Vector<amrex::MultiFab*> cost,
                                  std::string& knapsack_weight_type,
                                  int nbins, Real Rmax);

    ///
    /// Write out all particles in the format of the particle_input.dat files so we
    ///       can initialize a run with them
    ///
    void WriteAsciiFileForInit (const std::string& filename);

    ///
    /// Write out all particles at a level
    ///
    void writeAllAtLevel (int lev = 0);

    void split_particles (Real time);

  private:
    int nlev;
    int finest_level;
    bool reset_neighborhood;

    amrex::IntVect m_sorting_bin = amrex::IntVect::TheZeroVector();

    static int domain_bc[6];

    static void ReadStaticParameters ();

    bool p_verbose;
};
#endif
