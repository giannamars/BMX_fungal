//
//     Copyright (c) 2013 Battelle Memorial Institute
//     Licensed under modified BSD License. A copy of this license can be found
//     in the LICENSE file in the top level directory of this distribution.
//
#include <bmx_pc.H>
#include <bmx.H>
#include <bmx_chem_species_parms.H>

using namespace amrex;
using PinnedTile = typename BMXParticleContainer::ContainerLike<amrex::PinnedArenaAllocator>::ParticleTileType;

void BMXParticleContainer::WriteToAscii (const std::string& vtk_filename,int nstep, Real time)
{
  int nsegments, nparticles;
  CountParticleTypes(nparticles, nsegments);
  amrex::Print()<<"Nparticles: "<<nparticles<<" Nsegments: "<<nsegments<<std::endl;
  if (nsegments > 0) WriteSegmentsToAscii(vtk_filename, nstep, time, nsegments);
  if (nparticles > 0) WriteParticlesToAscii(vtk_filename, nstep, time, nparticles);
}

void BMXParticleContainer::WriteSegmentsToAscii (const std::string& vtk_filename,int nstep, Real time, const int nsegments)
{
  BL_PROFILE("BMXParticleContainer::WriteSegmentsToAscii");

  // Number of vertices per particle
  int NPC = 2;
   
  const std::string& filename = amrex::Concatenate(vtk_filename,nstep,6);

  Real l_max_vol  = SPECIES::max_vol;
  Real particle_rad = pow((3.0*l_max_vol/(4.0*M_PI)),1.0/3.0);

  amrex::Print()  << " "                                              << std::endl;
  amrex::Print()  << "  Writing vtkfile " << filename <<  " at time " << time << std::endl;
  amrex::Print()  << " "                                              << std::endl;

  long total_number_of_particles = 0;

  for (int lev = 0; lev <= finest_level; lev++)
  {
      total_number_of_particles += NumberOfParticlesAtLevel(lev);
  }

  long total_number_of_vertices  = NPC * nsegments;
    
  // First write all the header information just once (from the I/O procesor)
  if (ParallelDescriptor::IOProcessor())
  {
     //
     // Have I/O processor open file and write out particle metadata.
     //
     std::ofstream File;
     File.open(filename + ".vtk", std::ios::out | std::ios::trunc);
   
     if (!File.good())
        amrex::FileOpenFailed(filename+".vtk");

     // vtkAscii File contents
     File << "# vtk DataFile Version 3.0 "         << std::endl;
     File << "hexagon using polygon function vtk " << std::endl;
     File << "ASCII"                               << std::endl;
     File << "DATASET POLYDATA"                    << std::endl;
     File << " "                                   << std::endl;
     File << "POINTS" << " " << total_number_of_vertices << " " << "float" << std::endl;
     File.flush();
     File.close();
     if (!File.good())
         amrex::Abort("BMXarticleContainer::WriteVTKAscii():  problem writing file");
  }

  ParallelDescriptor::Barrier();
 
  const int MyProc = ParallelDescriptor::MyProc();

  for (int proc = 0; proc < ParallelDescriptor::NProcs(); proc++)
  {
      if (MyProc == proc)
      {
          //
          // Each CPU opens the file for appending and adds its particles.
          //
          VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

          std::ofstream File;

          File.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

          File.open(filename+".vtk", std::ios::out|std::ios::app);

          File.precision(8);

          if (!File.good())
              amrex::FileOpenFailed(filename+".vtk");

          for (int lev = 0; lev <= finest_level; lev++)
          {
            const auto& plevel = GetParticles(lev);
            for (const auto& kv : plevel)
            {
              PinnedTile pinned_tile;
              pinned_tile.define(NumRuntimeRealComps(),NumRuntimeIntComps());
              pinned_tile.resize(kv.second.numParticles());
              amrex::copyParticles(pinned_tile, kv.second, 0, 0,
                  kv.second.numParticles());
              const auto& particles = pinned_tile.GetArrayOfStructs();

              for (int i = 0; i < particles.numParticles(); ++i)
              {
                const Real* par = &particles[i].rdata(0); 

                
                int type = particles[i].idata(intIdx::cell_type);
                Real x1, y1, z1, x2, y2, z2;
                Real theta = par[realIdx::theta];
                Real phi = par[realIdx::phi];
                Real ct = cos(theta);
                Real st = sin(theta);
                Real cp = cos(phi);
                Real sp = sin(phi);
                Real nx = st*cp;
                Real ny = st*sp;
                Real nz = ct;
                if (type == cellType::FUNGI) {
                  Real c_length = par[realIdx::c_length];
                  x1 = particles[i].pos(0) + 0.5*nx*c_length;
                  y1 = particles[i].pos(1) + 0.5*ny*c_length;
                  z1 = particles[i].pos(2) + 0.5*nz*c_length;

                  x2 = particles[i].pos(0) - 0.5*nx*c_length;
                  y2 = particles[i].pos(1) - 0.5*ny*c_length;
                  z2 = particles[i].pos(2) - 0.5*nz*c_length;

                  File << x1 << " " << y1 << " " << z1 << std::endl;
                  File << x2 << " " << y2 << " " << z2 << std::endl;
                  File << " " << std::endl;
                }
              } // i
            } // loop over particles at level "lev"
          } // levels
          File.flush();
          File.close();
          if (!File.good())
            amrex::Abort("BMXarticleContainer::WriteVTKAscii():  problem writing file");
        } // MyProc
        ParallelDescriptor::Barrier();
    } // loop over proc

    // POLYGONS # of cells #, of vertices per cell   
    if (ParallelDescriptor::IOProcessor())
    {
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::ofstream File;
        File.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        File.open(filename+".vtk", std::ios::out|std::ios::app);
        if (!File.good())
            amrex::FileOpenFailed(filename+".vtk");

        File << "POLYGONS" << " " <<       nsegments << " " << 
                                     (2+1)*nsegments << std::endl;

        // Write polygon info
        for (long i = 0; i < nsegments; ++i)
        {
             File << NPC << " "; 
             for (int j = 0; j < NPC; j++)
                File << i*NPC+j << " "; 
             File     << std::endl;
        }

        File << " "                                                       << std::endl;
        File << "CELL_DATA"     << " "       << nsegments << std::endl;
        File << "POINT_DATA"    << " "       << total_number_of_vertices  << std::endl; 
        File << "SCALARS "      << "cell "   << "int"                         << std::endl;
        File << "LOOKUP_TABLE " << "default"                                  << std::endl;

        // Write lookup info - this will be used to color the vertices of one cell the same color
        for (long i = 0; i < nsegments; ++i)
        {
             for (int j = 0; j < NPC; j++)
                File << i << " "; 
             File     << std::endl;
        }
        File     << std::endl;

        File << "SCALARS "      << "conc_B "   << "float"                         << std::endl;
        File << "LOOKUP_TABLE " << "default"                                  << std::endl;

        File.flush();
        File.close();
        if (!File.good())
             amrex::Abort("BMXarticleContainer::WriteVTKAscii():  problem writing file");
     } // I/O Proc

  ParallelDescriptor::Barrier();
 
  for (int proc = 0; proc < ParallelDescriptor::NProcs(); proc++)
  {
      if (MyProc == proc)
      {
          //
          // Each CPU opens the file for appending and adds its particles.
          //
          VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

          std::ofstream File;

          File.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

          File.open(filename+".vtk", std::ios::out|std::ios::app);

          File.precision(8);

          if (!File.good())
              amrex::FileOpenFailed(filename+".vtk");

          for (int lev = 0; lev <= finest_level; lev++)
          {
            const auto& plevel = GetParticles(lev);
            for (const auto& kv : plevel)
            {
              PinnedTile pinned_tile;
              pinned_tile.define(NumRuntimeRealComps(),NumRuntimeIntComps());
              pinned_tile.resize(kv.second.numParticles());
              amrex::copyParticles(pinned_tile, kv.second, 0, 0,
                  kv.second.numParticles());
              const auto& particles = pinned_tile.GetArrayOfStructs();

              for (int i = 0; i < particles.numParticles(); ++i)
              {
                const Real* par = &particles[i].rdata(0); 

                if (particles[i].idata(intIdx::cell_type) == cellType::FUNGI) {
                  Real concB = par[realIdx::first_data+1];

                  File << concB << " " << concB << std::endl;
                }
              } // i
            } // loop over particles at level "lev"
          } // levels
          File.flush();
          File.close();
          if (!File.good())
             amrex::Abort("BMXarticleContainer::WriteVTKAscii():  problem writing file");
        } // MyProc
        ParallelDescriptor::Barrier();
    } // loop over proc

}

void BMXParticleContainer::WriteParticlesToAscii (const std::string& vtk_filename,int nstep, Real time, const int nsegments)
{
  BL_PROFILE("BMXParticleContainer::WriteSegmentsToAscii");

  // Number of vertices per particle
  int NPC = 2;
   
  std::string particlefile = vtk_filename;
  particlefile.append("_part");
  const std::string& filename = amrex::Concatenate(particlefile,nstep,6);

  Real l_max_vol  = SPECIES::max_vol;
  Real particle_rad = pow((3.0*l_max_vol/(4.0*M_PI)),1.0/3.0);

  amrex::Print()  << " "                                              << std::endl;
  amrex::Print()  << "  Writing vtkfile " << filename <<  " at time " << time << std::endl;
  amrex::Print()  << " "                                              << std::endl;

  long total_number_of_particles = 0;

  for (int lev = 0; lev <= finest_level; lev++)
  {
      total_number_of_particles += NumberOfParticlesAtLevel(lev);
  }

  long total_number_of_vertices  = NPC * nsegments;
    
  // First write all the header information just once (from the I/O procesor)
  if (ParallelDescriptor::IOProcessor())
  {
     //
     // Have I/O processor open file and write out particle metadata.
     //
     std::ofstream File;
     File.open(filename + ".vtk", std::ios::out | std::ios::trunc);
   
     if (!File.good())
        amrex::FileOpenFailed(filename+".vtk");

     // vtkAscii File contents
     File << "# vtk DataFile Version 3.0 "         << std::endl;
     File << "hexagon using polygon function vtk " << std::endl;
     File << "ASCII"                               << std::endl;
     File << "DATASET POLYDATA"                    << std::endl;
     File << " "                                   << std::endl;
     File << "POINTS" << " " << total_number_of_vertices << " " << "float" << std::endl;
     File.flush();
     File.close();
     if (!File.good())
         amrex::Abort("BMXarticleContainer::WriteVTKAscii():  problem writing file");
  }

  ParallelDescriptor::Barrier();
 
  const int MyProc = ParallelDescriptor::MyProc();

  for (int proc = 0; proc < ParallelDescriptor::NProcs(); proc++)
  {
      if (MyProc == proc)
      {
          //
          // Each CPU opens the file for appending and adds its particles.
          //
          VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

          std::ofstream File;

          File.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

          File.open(filename+".vtk", std::ios::out|std::ios::app);

          File.precision(8);

          if (!File.good())
              amrex::FileOpenFailed(filename+".vtk");

          for (int lev = 0; lev <= finest_level; lev++)
          {
            const auto& plevel = GetParticles(lev);
            for (const auto& kv : plevel)
            {
              PinnedTile pinned_tile;
              pinned_tile.define(NumRuntimeRealComps(),NumRuntimeIntComps());
              pinned_tile.resize(kv.second.numParticles());
              amrex::copyParticles(pinned_tile, kv.second, 0, 0,
                  kv.second.numParticles());
              const auto& particles = pinned_tile.GetArrayOfStructs();

              for (int i = 0; i < particles.numParticles(); ++i)
              {
                const Real* par = &particles[i].rdata(0); 

                
                int type = particles[i].idata(intIdx::cell_type);
                Real x1, y1, z1, x2, y2, z2;
                Real theta = par[realIdx::theta];
                Real phi = par[realIdx::phi];
                Real ct = cos(theta);
                Real st = sin(theta);
                Real cp = cos(phi);
                Real sp = sin(phi);
                Real nx = st*cp;
                Real ny = st*sp;
                Real nz = ct;
                if (type == cellType::YEAST) {
                  Real c_length = 0.01*par[realIdx::radius];
                  x1 = particles[i].pos(0) + 0.5*nx*c_length;
                  y1 = particles[i].pos(1) + 0.5*ny*c_length;
                  z1 = particles[i].pos(2) + 0.5*nz*c_length;

                  x2 = particles[i].pos(0) - 0.5*nx*c_length;
                  y2 = particles[i].pos(1) - 0.5*ny*c_length;
                  z2 = particles[i].pos(2) - 0.5*nz*c_length;

                  File << x1 << " " << y1 << " " << z1 << std::endl;
                  File << x2 << " " << y2 << " " << z2 << std::endl;
                  File << " " << std::endl;
                }
              } // i
            } // loop over particles at level "lev"
          } // levels
          File.flush();
          File.close();
          if (!File.good())
            amrex::Abort("BMXarticleContainer::WriteVTKAscii():  problem writing file");
        } // MyProc
        ParallelDescriptor::Barrier();
    } // loop over proc

    // POLYGONS # of cells #, of vertices per cell   
    if (ParallelDescriptor::IOProcessor())
    {
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::ofstream File;
        File.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        File.open(filename+".vtk", std::ios::out|std::ios::app);
        if (!File.good())
            amrex::FileOpenFailed(filename+".vtk");

        File << "POLYGONS" << " " <<       nsegments << " " << 
                                     (2+1)*nsegments << std::endl;

        // Write polygon info
        for (long i = 0; i < nsegments; ++i)
        {
             File << NPC << " "; 
             for (int j = 0; j < NPC; j++)
                File << i*NPC+j << " "; 
             File     << std::endl;
        }

        File << " "                                                       << std::endl;
        File << "CELL_DATA"     << " "       << nsegments << std::endl;
        File << "POINT_DATA"    << " "       << total_number_of_vertices  << std::endl; 
        File << "SCALARS "      << "cell "   << "int"                         << std::endl;
        File << "LOOKUP_TABLE " << "default"                                  << std::endl;

        // Write lookup info - this will be used to color the vertices of one cell the same color
        for (long i = 0; i < nsegments; ++i)
        {
             for (int j = 0; j < NPC; j++)
                File << i << " "; 
             File     << std::endl;
        }
        File     << std::endl;

        File << "SCALARS "      << "conc_E "   << "float"                         << std::endl;
        File << "LOOKUP_TABLE " << "default"                                  << std::endl;

        File.flush();
        File.close();
        if (!File.good())
             amrex::Abort("BMXarticleContainer::WriteVTKAscii():  problem writing file");
     } // I/O Proc

  ParallelDescriptor::Barrier();
 
  for (int proc = 0; proc < ParallelDescriptor::NProcs(); proc++)
  {
      if (MyProc == proc)
      {
          //
          // Each CPU opens the file for appending and adds its particles.
          //
          VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

          std::ofstream File;

          File.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

          File.open(filename+".vtk", std::ios::out|std::ios::app);

          File.precision(8);

          if (!File.good())
              amrex::FileOpenFailed(filename+".vtk");

          for (int lev = 0; lev <= finest_level; lev++)
          {
            const auto& plevel = GetParticles(lev);
            for (const auto& kv : plevel)
            {
              PinnedTile pinned_tile;
              pinned_tile.define(NumRuntimeRealComps(),NumRuntimeIntComps());
              pinned_tile.resize(kv.second.numParticles());
              amrex::copyParticles(pinned_tile, kv.second, 0, 0,
                  kv.second.numParticles());
              const auto& particles = pinned_tile.GetArrayOfStructs();

              for (int i = 0; i < particles.numParticles(); ++i)
              {
                const Real* par = &particles[i].rdata(0); 

                if (particles[i].idata(intIdx::cell_type) == cellType::YEAST) {
                  Real concE = par[realIdx::first_data+4];

                  File << concE << " " << concE << std::endl;
                }
              } // i
            } // loop over particles at level "lev"
          } // levels
          File.flush();
          File.close();
          if (!File.good())
             amrex::Abort("BMXarticleContainer::WriteVTKAscii():  problem writing file");
        } // MyProc
        ParallelDescriptor::Barrier();
    } // loop over proc

}

void BMXParticleContainer::CountParticleTypes (int &nparticles, int &nsegments)
{
  BL_PROFILE("BMXParticleContainer::WriteSegmentsToAscii");

#if 1
  amrex::ReduceOps<ReduceOpSum, ReduceOpSum> reduce_ops;
  auto r = amrex::ParticleReduce<ReduceData<int, int>> (
    *this, [=] AMREX_GPU_DEVICE (const ParticleType& p) noexcept -> amrex::GpuTuple<int, int>
    {
      int type = p.idata(intIdx::cell_type);
      if (type == cellType::FUNGI) {
      return {0, 1};
    } else {
      return {1, 0};
    }
    }, reduce_ops);

  nparticles = amrex::get<0>(r);
  nsegments = amrex::get<1>(r);
#else
  nparticles = 0;
  nsegments = 0;
    
  const int MyProc = ParallelDescriptor::MyProc();

  for (int lev = 0; lev <= finest_level; lev++)
  {
    int n_at_lev = this->NumberOfParticlesAtLevel(lev);

    if (n_at_lev == 0) continue;
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (BMXParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      PairIndex index(pti.index(), pti.LocalTileIndex());

      auto& plevel = GetParticles(lev);
      auto& ptile = plevel[index];
      auto& aos   = ptile.GetArrayOfStructs();
      ParticleType* pstruct = aos().dataPtr();


      Gpu::DeviceScalar<int> npart_gpu(0);
      Gpu::DeviceScalar<int> nseg_gpu(0);
      int* npart = npart_gpu.dataPtr();
      int* nseg = nseg_gpu.dataPtr();

      const int nrp = GetParticles(lev)[index].numRealParticles();

      amrex::ParallelFor(nrp, [pstruct,npart,nseg,MyProc]
          AMREX_GPU_DEVICE(int i) noexcept
      {
        auto& particle = pstruct[i];

        int type = particle.idata(intIdx::cell_type);
        if (type == cellType::FUNGI) {
          (*nseg)++;
          std::printf("p[%d] Increment nseg: %d\n",MyProc,*nseg);
        } else {
          (*npart)++;
        }
      }); // loop over particles at level "lev"
      nparticles += npart_gpu.dataValue();
      nsegments += nseg_gpu.dataValue();
    }
  } // levels
#endif
  int buf[2];
  buf[0] = nparticles;
  buf[1] = nsegments;
  ParallelDescriptor::ReduceIntSum(buf,2);
  nparticles = buf[0];
  nsegments = buf[1];
}
