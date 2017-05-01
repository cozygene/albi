import sys
if sys.version_info[0] == 3:
	from .albi_lib import estimate_distributions_eigenvectors
	from .albi_lib import estimate_distributions_general
	from .albi_lib import build_heritability_cis

	from .fiesta_lib import calculate_cis_eigenvectors
	from .fiesta_lib import calculate_cis_general
else:
	from albi_lib import estimate_distributions_eigenvectors
	from albi_lib import estimate_distributions_general
	from albi_lib import build_heritability_cis

	from fiesta_lib import calculate_cis_eigenvectors
	from fiesta_lib import calculate_cis_general


__doc__ = "ALBI"