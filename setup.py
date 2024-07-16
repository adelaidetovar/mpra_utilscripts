from setuptools import setup

setup(
    name = "mpra_utilscripts",
    packages = ['bc_count', 'subasm'],

    entry_points = {
	    'console_scripts': [ 'filter_starcode_bcgroups = bc_count.filter_starcode_bcgroups:main',
		                     'join_starcode_histos = bc_count.join_starcode_histos:main',
		                     'jointhisto_to_sampbysamp_mtx = bc_count.jointhisto_to_sampbysamp_mtx:main',
		                     'plot_tag_count_histos = bc_count.plot_tag_count_histos:main',
		                     'read_stats = bc_count.read_stats:main',
		                     'samp2samp_mtx_cluster_plots = bc_count.samp2samp_mtx_cluster_plots:main',
		                     'read_stats_subasm = subasm.read_stats_subasm:main' ],
		                     'bc_stats = bc_count.bc_stats:main'

	}
)