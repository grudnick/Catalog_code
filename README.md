# Catalog_code
Code I use to manipulate and match catalogs

catalog_match_trans.py is an obsolete file.  It has been replaced by catalog_match_sky_trans.py, which is set up to work on catalogs in sky coordinates only and catalog_match_im_trans.py, which is designed to work on image coordinates.

These are called by image_coord_transform_im.py which contains the master cat_im_match() routine.

To run the code you need to go into a python environment and:
import image_coord_transform_im as pcti
septol = <matching tolerance in pixels>
#this matches the input and reference catalog and runs geomap and
#geoxytran on these catalogs
pcti.cat_match_im(septol, im_trans = 1)
  
Read the header of image_coord_transform_im.py for more information
