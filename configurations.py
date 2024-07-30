import os
home_directory = os.path.expanduser("~")
CONFIG = {
    1: {
        "output_dir": 'HSX',
        "wout": 'wout_HSX.nc',
        "params": { 'nphi': 121,'nlambda': 25,'nperiod': 2.0,'nstep': 350,'dt': 0.4,
                    'aky_min': 0.3,'aky_max': 3.0,'naky': 8,'LN': 1.0,'LT': 3.0,
                    's_radius': 0.25,'alpha_fieldline': 0,'ngauss': 3,'negrid': 8,'vnewk': 0.01
                  },
    },
    2: {
        "output_dir": 'W7-X',
        "wout": 'wout_W7-X.nc',
        "params": { 'nphi': 89,'nlambda': 29,'nperiod': 1.2,'nstep': 350,'dt': 0.4,
                    'aky_min': 0.3,'aky_max': 3.0,'naky': 8,'LN': 1.0,'LT': 3.0,
                    's_radius': 0.25,'alpha_fieldline': 0,'ngauss': 3,'negrid': 8,'vnewk': 0.01
                  },
    },
    3: {
        "output_dir": 'nfp1_QI',
        "wout": 'wout_nfp1_QI.nc',
        "params": { 'nphi': 101,'nlambda': 23,'nperiod': 2.0,'nstep': 300,'dt': 0.4,
                    'aky_min': 0.3,'aky_max': 4.0,'naky': 8,'LN': 1.0,'LT': 3.0,
                    's_radius': 0.25,'alpha_fieldline': 0,'ngauss': 3,'negrid': 8,'vnewk': 0.01
                  },
    },
    4: {
        "output_dir": 'nfp4_QH',
        "wout": 'wout_nfp4_QH.nc',
        "params": { 'nphi': 121,'nlambda': 25,'nperiod': 2.5,'nstep': 350,'dt': 0.4,
                    'aky_min': 0.3,'aky_max': 3.0,'naky': 6,'LN': 1.0,'LT': 3.0,
                    's_radius': 0.25,'alpha_fieldline': 0,'ngauss': 3,'negrid': 8,'vnewk': 0.01
                  },
    },
    5: {
        "output_dir": 'nfp2_QA',
        "wout": 'wout_nfp2_QA.nc',
        "params": { 'nphi': 89,'nlambda': 25,'nperiod': 3.0,'nstep': 280,'dt': 0.4,
                    'aky_min': 0.4,'aky_max': 3.0,'naky': 6,'LN': 1.0,'LT': 3.0,
                    's_radius': 0.25,'alpha_fieldline': 0,'ngauss': 3,'negrid': 8,'vnewk': 0.01
                  },
    }
}