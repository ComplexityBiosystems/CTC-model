# python executables
MAKE_BODYSYSTEM=code/make_bodysystem.py
FIND_CONTACT_POINTS=code/find_contact_points.py
STITCH_ARTERIES_VEINS=code/stitch_arteries_veins.py
LAUNCH_TRACERS=code/launch_tracers_fullsystem.py
SOLVE_FLOWS_NETWORK=code/solve_flows_network.py
OPTIMIZE_ORGAN_COEFFICIENTS=code/optimize_organ_coefficients.py

# artery data
ARTERY_PICKLED=data/SelectedBodyParts3D/artery/pickled
ARTERY_PLUS_BRIDGES=data/SelectedBodyParts3D/artery/artery_plus_bodybridges.txt
ARTERY_MINUS_BRIDGES=data/SelectedBodyParts3D/artery/artery_minus_bodybridges.txt
ARTERY_BODYSYSTEM=output/data/bodysystems/artery/artery_bodysystem.p

# vein data
VEIN_PICKLED=data/SelectedBodyParts3D/vein/pickled
VEIN_PLUS_BRIDGES=data/SelectedBodyParts3D/vein/vein_plus_bodybridges.txt
VEIN_MINUS_BRIDGES=data/SelectedBodyParts3D/vein/vein_minus_bodybridges.txt
VEIN_BODYSYSTEM=output/data/bodysystems/vein/vein_bodysystem.p

# flow data
FLOW_NETWORKS_DIR=output/data/flow_networks
OPTIMIZED_MAIN_NETWORK=output/data/flow_networks/optimized_main_network.p
SOLVED_MAIN_NETWORK=output/data/solved_flow_networks/solved_main_network.p
SOLVED_PULMONARY_NETWORK=output/data/solved_flow_networks/solved_pulmonary_network.p
SOLVED_TOTAL_CLOSED_NETWORK=output/data/solved_flow_networks/solved_total_closed.p


# boundary conditions
MAIN_BOUNDARY_CONDITIONS=output/data/boundary_conditions/main_boundary_conditions.csv
PULMONARY_BOUNDARY_CONDITIONS=output/data/boundary_conditions/pulmonary_boundary_conditions.csv

# artery-vein connections
MAIN_CONNECTIONS=output/data/flow_networks/main_connections.p
PULMONARY_CONNECTIONS=output/data/flow_networks/main_connections.p

# organ data
ORGAN_OBJ_FILES=data/SelectedBodyParts3D/organ/*/*.obj
ORGAN_LIST=brain large_intestine left_lung pancreas right_lung stomach heart left_kidney liver right_kidney small_intestine prostate gallbladder urinary_bladder left_adrenal_gland right_adrenal_gland

# contact points
CONTACT_POINTS_FILES=data/SelectedBodyParts3D/organ/*/*_*_bodysystem_solved_contacts.txt

# aliases
solve_main_system: ${SOLVED_MAIN_NETWORK}
solve_pulmonary_system: ${SOLVED_PULMONARY_NETWORK}


# optimize organ coefficients
${OPTIMIZED_MAIN_NETWORK}: ${OPTIMIZE_ORGAN_COEFFICIENTS} ${FLOW_NETWORKS_DIR}/main_network.p ${MAIN_BOUNDARY_CONDITIONS} ${FLOW_NETWORKS_DIR}/main_connections.p output/data/optimization/organs_dict.p output/data/optimization/organs_targets.p
	python ${OPTIMIZE_ORGAN_COEFFICIENTS} \
	--input_network ${FLOW_NETWORKS_DIR}/main_network.p \
	--output_network ${OPTIMIZED_MAIN_NETWORK} \
	--output_coefs ${FLOW_NETWORKS_DIR}/optimized_main_network_coefs.txt \
	--boundary_conditions ${MAIN_BOUNDARY_CONDITIONS} \
	--connections ${FLOW_NETWORKS_DIR}/main_connections.p \
	--organs_nodes_dict output/data/optimization/organs_dict.p \
	--organs_flow_fraction_target output/data/optimization/organs_targets.p \
	-v \
	--maxiter 50

# solve flows main system
${SOLVED_MAIN_NETWORK}: ${SOLVE_FLOWS_NETWORK} ${OPTIMIZED_MAIN_NETWORK} ${MAIN_BOUNDARY_CONDITIONS}
	python ${SOLVE_FLOWS_NETWORK} \
	--network_pickled ${OPTIMIZED_MAIN_NETWORK} \
	--boundary_conditions_csv ${MAIN_BOUNDARY_CONDITIONS} \
	--output_file ${SOLVED_MAIN_NETWORK}

# solve flows pulmonary system
${SOLVED_PULMONARY_NETWORK}: ${SOLVE_FLOWS_NETWORK} ${FLOW_NETWORKS_DIR}/pulmonary_* ${PULMONARY_BOUNDARY_CONDITIONS}
	python ${SOLVE_FLOWS_NETWORK} \
	--network_pickled ${FLOW_NETWORKS_DIR}/pulmonary_network.p \
	--boundary_conditions_csv ${PULMONARY_BOUNDARY_CONDITIONS} \
	--output_file ${SOLVED_PULMONARY_NETWORK}

# find contact points
contact_points: ${FIND_CONTACT_POINTS} ${ORGAN_OBJ_FILES}
	for ORGAN in ${ORGAN_LIST} ; do \
		python ${FIND_CONTACT_POINTS} \
		--solved_flows_network ${SOLVED_MAIN_NETWORK} \
		--organ data/SelectedBodyParts3D/organ/$$ORGAN/$$ORGAN.obj \
		--output_dir output/data/organs_contact_points \
		--tolerance 10 \
		--verbose ; \
		echo ; \
	done

# launch CTC simulations (fast, not storing full trajectory)
launch_tracers: ${LAUNCH_TRACERS} ${SOLVED_TOTAL_CLOSED_NETWORK} ${ORGAN_OBJ_FILES}
	for SSIZE in 1000 ; do \
		for SUFFIX in $$(seq -w 0 0) ; do \
			for ORGAN in ${ORGAN_LIST} ; do \
				python ${LAUNCH_TRACERS} \
				--drop_trajectories \
				--organ $$ORGAN \
				--suffix _notrajs_chunk$$SUFFIX \
				--network ${SOLVED_TOTAL_CLOSED_NETWORK} \
				-n $$SSIZE \
				--output_dir output/data/tracers_notrajs/ & \
				echo & \
			done ; \
		done ; \
	done

# launch CTC simulations (slower, storing everything)
launch_tracers_with_times: ${LAUNCH_TRACERS} ${SOLVED_TOTAL_CLOSED_NETWORK} ${ORGAN_OBJ_FILES}
	for SSIZE in 500 ; do \
		for SUFFIX in $$(seq -w 0 0) ; do \
			for ORGAN in ${ORGAN_LIST} ; do \
				python ${LAUNCH_TRACERS} \
				--keep_times \
				--organ $$ORGAN \
				--suffix _withtimes_chunk$$SUFFIX \
				--network ${SOLVED_TOTAL_CLOSED_NETWORK} \
				-n $$SSIZE \
				--output_dir output/data/tracers_with_times/ & \
				echo & \
			done ; \
		done ; \
	done

