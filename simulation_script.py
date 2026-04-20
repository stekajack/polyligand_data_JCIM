from pressomancy.simulation import Simulation, Crowder, Filament, Quartet, Quadriplex
from pressomancy.helper_functions import BondWrapper
import espressomd
from espressomd.io.writer import vtf
from espressomd import checkpointing
import argparse
import os
import sys as sysos
import numpy as np
import time
import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)
end_of_time = time.time()+84600  # 23.5h in seconds
N_avog = 6.02214076e23
parser = argparse.ArgumentParser()
parser.add_argument('-no_obj', '--no_obj', type=int, required=True,
                    help='int number of quartets')
parser.add_argument('-no_crowders', '--no_crowders', type=int, required=True,
                    help='int number or crowders')
parser.add_argument('-no_per', '--no_per', type=int, required=True,
                    help='int N of particles per quartet')
parser.add_argument('-concentration', '--concentration', type=float, required=True,
                    help='concentration of particles in mM')
parser.add_argument('-vdW', '--vdW', type=float,
                    required=True, help='float vdW epsilon value')
parser.add_argument('-vdW_ligand', '--vdW_ligand', type=float,
                    required=True, help='float vdW epsilon value')
parser.add_argument('-part_per_filament', '--part_per_filament', type=int,
                    required=True, help='number of monomer per filament. If <=1, skips the filament contruction')
parser.add_argument('-part_per_ligand', '--part_per_ligand', type=int,
                    required=True, help='number of monomer per ligand. If <=1, skips the filament contruction')
parser.add_argument('-path_data', '--path_data', type=str, required=True,
                    help='absolute path to data')
parser.add_argument('-MODE', '--MODE', type=str, required=True,
                    help='do i start clean or load a checkpoint')
parser.add_argument('-bonding_mode', '--bonding_mode', type=str, required=True,
                    help='how do i bond G4s')
args = parser.parse_args()

context_string = f'_{args.no_obj}_{args.no_per}_{args.no_crowders}_{args.concentration}_{args.vdW}_{args.vdW_ligand}_{args.part_per_filament}_{args.part_per_ligand}'
logging.info(context_string)
logging.info(f'args.path_data: {args.path_data}')

dir_path = args.path_data
traj_path=os.path.join(dir_path, 'trajectory'+context_string+'.vtf')
h5_data_path=os.path.join(dir_path, 'custom_data_wip'+context_string+'.h5')
logging.info(traj_path)

sigma = 1.
rho_si = args.concentration*N_avog
N = int(args.no_obj/3)
vol = N/rho_si
box_l = pow(vol, 1/3)
_box_l = box_l/0.4e-09
box_dim = _box_l*np.ones(3)
_rho = N/pow(_box_l, 3)

sheets_per_quad = 3

logging.info(f'box_dim: {box_dim}', )
sim_inst = Simulation(box_dim=box_dim)
sim_inst.set_sys(timestep=0.005, min_global_cut=7.)

if args.MODE == 'LOAD_NEW':
    # //////////////////////////////////////////////////////////////////////////////////////////////
    print(f'Loading checkpoint checkpoint {context_string} at {dir_path}')
    checkpoint = checkpointing.Checkpoint(
        checkpoint_id=f"checkpoints_{context_string}", checkpoint_path=dir_path)
    try:
        checkpoint.load()
    except Exception:
        exc_type, value, traceback = sysos.exc_info()
        logging.info("Failed with exception [%s,%s ,%s]" %
              (exc_type, value, traceback))
        sim_inst.sys.part.clear()
        logging.info('Retrying to load')
        checkpoint.load()

    fp = open(traj_path, mode='a')

if args.MODE == 'NEW':

    quartet_configuration = Quartet.config.specify(sigma=sigma, n_parts=args.no_per, type='solid', espresso_handle=sim_inst.sys)
    quartets = [Quartet(config=quartet_configuration) for x in range(args.no_obj)]
    sim_inst.store_objects(quartets)

    bond_quad = BondWrapper(espressomd.interactions.FeneBond(k=10., r_0=2., d_r_max=2*1.5))
    grouped_quartets = [quartets[i:i+sheets_per_quad]
                        for i in range(0, len(quartets), sheets_per_quad)]
    # size estimated as the radius of circumscribed sphere!
    quadriplex_configuration_list = [Quadriplex.config.specify(
        espresso_handle=sim_inst.sys,associated_objects=elem,bonding_mode=args.bonding_mode,size=np.sqrt(3)*5.13,bond_handle=bond_quad) for elem in grouped_quartets]
    quadriplex = [Quadriplex(config=elem) for elem in quadriplex_configuration_list]
    sim_inst.store_objects(quadriplex)

    if args.part_per_filament > 1:
        bond_filam = BondWrapper(espressomd.interactions.FeneBond(k=10., r_0=2., d_r_max=2*1.5))
        grouped_quadriplexes = [quadriplex[i:i+args.part_per_filament:]
                                for i in range(0, len(quadriplex), args.part_per_filament)]
        a=5.13
        b=args.part_per_filament*4.+bond_filam.r_0*(args.part_per_filament-1)+1.13
        d=np.sqrt(pow(np.sqrt(2)*a,2)+b**2)
        filament_configuration_list = [Filament.config.specify(size=d, n_parts=args.part_per_filament, espresso_handle=sim_inst.sys, bond_handle=bond_filam, associated_objects=elem, spacing=6.) for elem in grouped_quadriplexes]
        filaments = [Filament(config=configuration) for configuration in filament_configuration_list]
        sim_inst.store_objects(filaments)
        sim_inst.set_objects(filaments)
        for filament in filaments:
            filament.bond_quadriplexes()
    else:
        sim_inst.set_objects(quadriplex)

    if args.vdW:
        for el in quadriplex:
            el.add_patches_triples()
        sim_inst.set_vdW(key=('patch',), lj_eps=args.vdW, lj_size=2.)

    if args.no_crowders:
        sim_inst.sys.integrator.run(0)
        crowder_configuration = Crowder.config.specify(size=pow(2,1/6), espresso_handle=sim_inst.sys)
        crowders = [Crowder(config=crowder_configuration)
                    for x in range(args.no_crowders)]
        sim_inst.store_objects(crowders)

        if args.part_per_ligand > 1:
            grouped_crowders = [crowders[i:i+args.part_per_ligand]
                                for i in range(0, len(crowders), args.part_per_ligand)]
            bender_pass = BondWrapper(espressomd.interactions.FeneBond(
            k=10, r_0=6, d_r_max=6*1.5))
            filament_configuration_list = [Filament.config.specify(
                n_parts=args.part_per_ligand, espresso_handle=sim_inst.sys, associated_objects=elem, size=6*args.part_per_ligand,bond_handle=bender_pass) for elem in grouped_crowders]
            filaments = [Filament(config=elem) for elem in filament_configuration_list]
            sim_inst.store_objects(filaments)
            sim_inst.set_objects(filaments)
            angle_harmonic=espressomd.interactions.AngleHarmonic(bend=10., phi0=np.pi)
            for filament in filaments:
                filament.bond_center_to_center(type_key='crowder')
                # filament.add_bending_potential(type_key='crowder',bond_handle=angle_harmonic)
        else:
            sim_inst.set_objects(crowders)

        sim_inst.set_steric(key=('real', 'virt','crowder'), wca_eps=1.)
        
        if args.vdW_ligand:
            sim_inst.set_vdW_custom(pairs=[('patch','crowder'),], lj_eps=[args.vdW_ligand,], lj_size=[1.,])
    else:
        sim_inst.set_steric(key=('real', 'virt'), wca_eps=1.)

    sim_inst.sys.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=sim_inst.seed)

    fp = open(traj_path, mode='w+t')
    vtf.writevsf(sim_inst.sys, fp)
    vtf.writevcf(sim_inst.sys, fp)
    fp.flush()

    # sim_inst.avoid_explosion(F_TOL=1e-2)
    # //////////////////////////////////////////////////////////////////////////////////////////////
    checkpoint = checkpointing.Checkpoint(
        checkpoint_id=f"checkpoints_{context_string}", checkpoint_path=dir_path)
    checkpoint.register("sim_inst.sys")
# /////////////////////////////////////////////////////////////////////////////////////////////


    # EQUILIBRATION_ITERATIONS = 5
    # EQUILIBRATION_INTERVAL = int(5e05)
    # print(sim_inst.sys.analysis.energy())
    # for el in range(EQUILIBRATION_ITERATIONS):
    #     logging.info(f'EQUILIBRATION_ITERATIONS: {el}')
    #     sim_inst.sys.integrator.run(EQUILIBRATION_INTERVAL)
    #     vtf.writevcf(sim_inst.sys, fp)

GLOBAL_COUNTER=sim_inst.inscribe_part_group_to_h5(group_type=[Filament,Crowder], h5_data_path=h5_data_path, mode=args.MODE)
benchmark_SAMPLING_INTERVAL = [0.,]
t1=0.
SAMPLING_ITERATIONS = 15
SAMPLING_INTERVAL = int(5e05)
while GLOBAL_COUNTER <= SAMPLING_ITERATIONS:
    t0 = time.time()
    if end_of_time-t1 > benchmark_SAMPLING_INTERVAL[-1]*2:
        logging.info(f"SAMPLING_ITERATIONS: {GLOBAL_COUNTER}" )
        sim_inst.sys.integrator.run(SAMPLING_INTERVAL)
        vtf.writevcf(sim_inst.sys, fp)
        fp.flush()
        sim_inst.write_part_group_to_h5(time_step=GLOBAL_COUNTER)

        if GLOBAL_COUNTER == SAMPLING_ITERATIONS:
            checkpoint.save()
            logging.info("final checkpoint.save() ran sucessfully. Exiting simulation")
        GLOBAL_COUNTER += 1

    else:
        checkpoint.save()
        logging.info("checkpoint.save() ran sucessfully. Exiting loop: ")
        break
    t1 = time.time()
    benchmark_SAMPLING_INTERVAL.append(t1-t0)

logging.info(f"benchmark SAMPLING_ITERATIONS: {np.mean(benchmark_SAMPLING_INTERVAL[1:])}")
