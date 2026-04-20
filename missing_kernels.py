import numpy as np
from itertools import pairwise
from pressomancy.analysis import H5DataSelector
import h5py
from pmtools.refractored_toolbox import determine_key_val_from_filename, check_breakage, unbreak_graph

def Ree_h5(data_path, template_hndl, box_dim, chunk=(-5, None, 1), norm=1., crit=1.47, extra_flag=None):
    data_with_context = {}
    data_file = h5py.File(data_path, "r")
    data = H5DataSelector(data_file, particle_group="Filament")
    monomer_no = int(determine_key_val_from_filename(template_hndl, data_path, 'what_monomer_number'))

    accumulated_ree = []
    start, end, step = chunk

    for col in data.timestep[start:end:step].timestep:
        filtered_fil_ids = []
        for myed in list(col.get_connectivity_values('Filament')):
            parts = col.select_particles_by_object('Filament', myed)
            types = parts.type.flatten()
            if 5 not in types:
                filtered_fil_ids.append(myed)

        filtered_fil_ids.sort()
        pf_indices = []
        filtered_pos = []

        for myed in filtered_fil_ids:
            ids_shuffled = col.select_particles_by_object('Filament', myed).id.flatten()
            pos_shuffled = col.select_particles_by_object('Filament', myed).pos

            order = np.argsort(ids_shuffled)
            ids_ordered = ids_shuffled[order]
            pos_ordered = pos_shuffled[order]

            # Remove patches so the remaining coordinates are only the filament beads.
            pf_indices.append(ids_ordered[:-2 * monomer_no])
            filtered_pos.extend(pos_ordered[:-2 * monomer_no])

        edges = [(int(x), int(y)) for pf_el in pf_indices for x, y in pairwise(pf_el)]
        g2 = ig.Graph(n=len(filtered_pos), edges=edges)
        g2.vs["pos"] = filtered_pos
        g2.simplify()
        decomposition = g2.decompose()

        for subgraph in decomposition:
            flag, pass_graph = check_breakage(subgraph, box_dim)
            if not flag:
                positions = unbreak_graph(pass_graph, box_dim)
            else:
                positions = np.array(subgraph.vs["pos"])

            com_pos = np.mean(positions.reshape(monomer_no, -1, 3), axis=1)
            ree = np.linalg.norm(com_pos[-1] - com_pos[0])
            accumulated_ree.append(ree)

    data_with_context[data_path] = accumulated_ree
    return data_with_context

def segments_h5(data_path, template_hndl, box_dim, chunk=(-5, None, 1), norm=1., crit=1.47, extra_flag=None):

    data_with_context = {}
    data_file=h5py.File(data_path, "r")
    data=H5DataSelector(data_file,particle_group="Filament")
    monomer_no = int(determine_key_val_from_filename(template_hndl,data_path,'what_monomer_number'))
    accumulated_lp_seg = []
    start, end, step = chunk
    for col in data.timestep[start:end:step].timestep:
        fitered_fil_ids=[]
        for myed in list(col.get_connectivity_values('Filament')):
            parts=col.select_particles_by_object('Filament',myed)
            types=parts.type.flatten()
            if 5 not in types:
                fitered_fil_ids.append(myed)

        fitered_fil_ids.sort()
        pf_indices=[]
        filtered_pos=[]
        for myed in fitered_fil_ids:
            ids_shuffled=col.select_particles_by_object('Filament',myed).id.flatten()
            pos_shuffled=col.select_particles_by_object('Filament',myed).pos
            order = np.argsort(ids_shuffled)
            ids_ordered = ids_shuffled[order]
            pos_ordered = pos_shuffled[order]
            # Remove patches for the iGraph unfolding to work correctly. 
            # Indices must be increasing monotonically (patches do not w.r.t rest of part)
            pf_indices.append(ids_ordered[:-2*monomer_no])
            filtered_pos.extend(pos_ordered[:-2*monomer_no])

        edges = [(int(x), int(y)) for pf_el in pf_indices for x,
                     y in pairwise(pf_el)]
        g2 = ig.Graph(n=len(filtered_pos), edges=edges)
        g2.vs["pos"] = filtered_pos
        g2.simplify()
        decomposition = g2.decompose()
        for subgraph in decomposition:
            flag, pass_graph = check_breakage(
                subgraph, box_dim)
            if not flag:
                positions = unbreak_graph(
                    pass_graph, box_dim)
            else:
                positions = np.array(subgraph.vs['pos'])

            com_pos = np.mean(positions.reshape(
                monomer_no, -1, 3), axis=1)
            ete_vec = com_pos[-1]-com_pos[0]
            segments = np.diff(com_pos, axis=0)
            seg_norms = np.linalg.norm(segments, axis=1)
            accumulated_lp_seg.append(seg_norms)            
    data_with_context[data_path] = accumulated_lp_seg
    return data_with_context