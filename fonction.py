def smooth(filename, fwhm) :
	import healpy as hp
	import numpy as np
	map = hp.read_map(filename)
	# map_size = hp.get_map_size(map)
	#hp.mollview(map, norm='hist')
	# print(hp.npix2nside(map_size))

        print fwhm
	fwhm_rad = np.radians(fwhm/60)
	print fwhm_rad
	map_gauss = hp.smoothing(map, fwhm_rad)
	#hp.mollview(map_gauss, norm='hist')
	print 'end of smooth'
	return map_gauss
