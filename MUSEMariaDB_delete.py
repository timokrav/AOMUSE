from pony.orm import *

# Create a database object from Pony
db = Database()

# The classes inherit db.Entity from Pony
class Target(db.Entity):
    #   ----- Attributes -----

    target_name = Required(str, unique=True)  # Required: Cannot be None

    #   ----- Relations -----

    exposures = Set('Exposure')  # One target contains a set of exposures

# Exposure table class
class Exposure(db.Entity):
    #   ----- Attributes -----

    observation_time = Required(str, unique=True)
    insMode = Required(str)
    datacube_header = Optional(LongStr)
    raw_exposure_header = Optional(LongStr)
    raw_exposure_data = Optional(LongStr)
    raw_exposure_filename = Optional(str, unique=True)
    raman_image_header = Optional(LongStr)
    pampelmuse_catalog = Optional(LongStr)
    pampelmuse_params = Optional(LongStr)
    prm_filename = Optional(str, unique=True)
    psf_table = Optional(LongStr)
    sources = Optional(LongStr)

    #   ----- Sky parameters -----
    sky_condition_start_time = Optional(float)
    sky_condition_start = Optional(LongStr)
    sky_comment_start = Optional(LongStr)
    sky_condition_end_time = Optional(float)
    sky_condition_end = Optional(LongStr)
    sky_comment_end = Optional(LongStr)

    #   ----- Relations -----

    target = Required('Target')  # One exposure belongs to a target


#   ----- Main -----

db.bind(provider='mysql', host='127.0.0.1', user='user', passwd='pass', db='dbname') # Establish the conection with the database
db.generate_mapping() # Map the classes with the database tables 
db.drop_all_tables(with_all_data = True) # Delete all the tables

