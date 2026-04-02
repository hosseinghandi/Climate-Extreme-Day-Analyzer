import math
import csv
import datetime
import ghpythonlib.components as gh

# Output containers
dry_bulb_temperature = []
relative_humidity = []
wind_speed = []
model_year = []



def epw_open(epw_path):
    if not epw_path:
        return None
    with open(epw_path, 'r') as epw:
        reader = csv.reader(epw)

        # --- LOCATION LINE ---
        location = next(reader)
        latitude = float(location[6])
        longitude = float(location[7])
        timezone = float(location[8])

        # --- SKIP HEADER LINES ---
        for _ in range(7):
            next(reader)

        # --- HOURLY DATA ---
        for row in reader:
            model_year.append(int(row[0]))
            dry_bulb_temperature.append(float(row[6]))
            relative_humidity.append(float(row[8]))
            wind_speed.append(float(row[21]))
    return (
        dry_bulb_temperature,
        wind_speed,
        relative_humidity
    )

def universal_thermal_climate_index(ta, vel, rh):
    # set upper and lower limits of air velocity according to Fiala model scenarios
    vel = 0.5 if vel < 0.5 else vel
    vel = 17 if vel > 17 else vel
    # metrics derived from the inputs used in the polynomial equation
    eh_pa = saturated_vapor_pressure_hpa(ta) * (rh / 100.0)  # partial vapor pressure
    pa_pr = eh_pa / 10.0  # convert vapour pressure to kPa
    d_tr = 0  # difference between radiant and air temperature
    # pre-calculate powers so we can re-use them
    utci_approx = _utci_polynomial(ta, d_tr, vel, pa_pr)
    return utci_approx

def saturated_vapor_pressure_hpa(db_temp):
    g = (-2836.5744, -6028.076559, 19.54263612, -0.02737830188, 0.000016261698,
         7.0229056e-10, -1.8680009e-13)
    tk = db_temp + 273.15  # air temp in K
    es = 2.7150305 * math.log(tk)
    for i, x in enumerate(g):
        es = es + (x * (tk**(i - 2)))
    es = math.exp(es) * 0.01
    return es

def _utci_polynomial(ta, d_tr, vel, pa_pr):
    # pre-calculate powers so we can re-use them
    ta2 = ta ** 2
    ta3 = ta ** 3
    ta4 = ta ** 4
    ta5 = ta ** 5
    ta6 = ta ** 6
    vel2 = vel ** 2
    vel3 = vel ** 3
    vel4 = vel ** 4
    vel5 = vel ** 5
    vel6 = vel ** 6
    d_tr2 = d_tr ** 2
    d_tr3 = d_tr ** 3
    d_tr4 = d_tr ** 4
    d_tr5 = d_tr ** 5
    d_tr6 = d_tr ** 6
    pa_pr2 = pa_pr ** 2
    pa_pr3 = pa_pr ** 3
    pa_pr4 = pa_pr ** 4
    pa_pr5 = pa_pr ** 5
    pa_pr6 = pa_pr ** 6
    # UTCI approximation calculation
    utci_approx = ta + \
        0.607562052 + \
        -0.0227712343 * ta + \
        8.06470249e-4 * ta2 + \
        -1.54271372e-4 * ta3 + \
        -3.24651735e-6 * ta4 + \
        7.32602852e-8 * ta5 + \
        1.35959073e-9 * ta6 + \
        -2.25836520 * vel + \
        0.0880326035 * ta * vel + \
        0.00216844454 * ta2 * vel + \
        -1.53347087e-5 * ta3 * vel + \
        -5.72983704e-7 * ta4 * vel + \
        -2.55090145e-9 * ta5 * vel + \
        -0.751269505 * vel2 + \
        -0.00408350271 * ta * vel2 + \
        -5.21670675e-5 * ta2 * vel2 + \
        1.94544667e-6 * ta3 * vel2 + \
        1.14099531e-8 * ta4 * vel2 + \
        0.158137256 * vel3 + \
        -6.57263143e-5 * ta * vel3 + \
        2.22697524e-7 * ta2 * vel3 + \
        -4.16117031e-8 * ta3 * vel3 + \
        -0.0127762753 * vel4 + \
        9.66891875e-6 * ta * vel4 + \
        2.52785852e-9 * ta2 * vel4 + \
        4.56306672e-4 * vel5 + \
        -1.74202546e-7 * ta * vel5 + \
        -5.91491269e-6 * vel6 + \
        0.398374029 * d_tr + \
        1.83945314e-4 * ta * d_tr + \
        -1.73754510e-4 * ta2 * d_tr + \
        -7.60781159e-7 * ta3 * d_tr + \
        3.77830287e-8 * ta4 * d_tr + \
        5.43079673e-10 * ta5 * d_tr + \
        -0.0200518269 * vel * d_tr + \
        8.92859837e-4 * ta * vel * d_tr + \
        3.45433048e-6 * ta2 * vel * d_tr + \
        -3.77925774e-7 * ta3 * vel * d_tr + \
        -1.69699377e-9 * ta4 * vel * d_tr + \
        1.69992415e-4 * vel2 * d_tr + \
        -4.99204314e-5 * ta * vel2 * d_tr + \
        2.47417178e-7 * ta2 * vel2 * d_tr + \
        1.07596466e-8 * ta3 * vel2 * d_tr + \
        8.49242932e-5 * vel3 * d_tr + \
        1.35191328e-6 * ta * vel3 * d_tr + \
        -6.21531254e-9 * ta2 * vel3 * d_tr + \
        -4.99410301e-6 * vel4 * d_tr + \
        -1.89489258e-8 * ta * vel4 * d_tr + \
        8.15300114e-8 * vel5 * d_tr + \
        7.55043090e-4 * d_tr2 + \
        -5.65095215e-5 * ta * d_tr2 + \
        -4.52166564e-7 * ta2 * d_tr2 + \
        2.46688878e-8 * ta3 * d_tr2 + \
        2.42674348e-10 * ta4 * d_tr2 + \
        1.54547250e-4 * vel * d_tr2 + \
        5.24110970e-6 * ta * vel * d_tr2 + \
        -8.75874982e-8 * ta2 * vel * d_tr2 + \
        -1.50743064e-9 * ta3 * vel * d_tr2 + \
        -1.56236307e-5 * vel2 * d_tr2 + \
        -1.33895614e-7 * ta * vel2 * d_tr2 + \
        2.49709824e-9 * ta2 * vel2 * d_tr2 + \
        6.51711721e-7 * vel3 * d_tr2 + \
        1.94960053e-9 * ta * vel3 * d_tr2 + \
        -1.00361113e-8 * vel4 * d_tr2 + \
        -1.21206673e-5 * d_tr3 + \
        -2.18203660e-7 * ta * d_tr3 + \
        7.51269482e-9 * ta2 * d_tr3 + \
        9.79063848e-11 * ta3 * d_tr3 + \
        1.25006734e-6 * vel * d_tr3 + \
        -1.81584736e-9 * ta * vel * d_tr3 + \
        -3.52197671e-10 * ta2 * vel * d_tr3 + \
        -3.36514630e-8 * vel2 * d_tr3 + \
        1.35908359e-10 * ta * vel2 * d_tr3 + \
        4.17032620e-10 * vel3 * d_tr3 + \
        -1.30369025e-9 * d_tr4 + \
        4.13908461e-10 * ta * d_tr4 + \
        9.22652254e-12 * ta2 * d_tr4 + \
        -5.08220384e-9 * vel * d_tr4 + \
        -2.24730961e-11 * ta * vel * d_tr4 + \
        1.17139133e-10 * vel2 * d_tr4 + \
        6.62154879e-10 * d_tr5 + \
        4.03863260e-13 * ta * d_tr5 + \
        1.95087203e-12 * vel * d_tr5 + \
        -4.73602469e-12 * d_tr6 + \
        5.12733497 * pa_pr + \
        -0.312788561 * ta * pa_pr + \
        -0.0196701861 * ta2 * pa_pr + \
        9.99690870e-4 * ta3 * pa_pr + \
        9.51738512e-6 * ta4 * pa_pr + \
        -4.66426341e-7 * ta5 * pa_pr + \
        0.548050612 * vel * pa_pr + \
        -0.00330552823 * ta * vel * pa_pr + \
        -0.00164119440 * ta2 * vel * pa_pr + \
        -5.16670694e-6 * ta3 * vel * pa_pr + \
        9.52692432e-7 * ta4 * vel * pa_pr + \
        -0.0429223622 * vel2 * pa_pr + \
        0.00500845667 * ta * vel2 * pa_pr + \
        1.00601257e-6 * ta2 * vel2 * pa_pr + \
        -1.81748644e-6 * ta3 * vel2 * pa_pr + \
        -1.25813502e-3 * vel3 * pa_pr + \
        -1.79330391e-4 * ta * vel3 * pa_pr + \
        2.34994441e-6 * ta2 * vel3 * pa_pr + \
        1.29735808e-4 * vel4 * pa_pr + \
        1.29064870e-6 * ta * vel4 * pa_pr + \
        -2.28558686e-6 * vel5 * pa_pr + \
        -0.0369476348 * d_tr * pa_pr + \
        0.00162325322 * ta * d_tr * pa_pr + \
        -3.14279680e-5 * ta2 * d_tr * pa_pr + \
        2.59835559e-6 * ta3 * d_tr * pa_pr + \
        -4.77136523e-8 * ta4 * d_tr * pa_pr + \
        8.64203390e-3 * vel * d_tr * pa_pr + \
        -6.87405181e-4 * ta * vel * d_tr * pa_pr + \
        -9.13863872e-6 * ta2 * vel * d_tr * pa_pr + \
        5.15916806e-7 * ta3 * vel * d_tr * pa_pr + \
        -3.59217476e-5 * vel2 * d_tr * pa_pr + \
        3.28696511e-5 * ta * vel2 * d_tr * pa_pr + \
        -7.10542454e-7 * ta2 * vel2 * d_tr * pa_pr + \
        -1.24382300e-5 * vel3 * d_tr * pa_pr + \
        -7.38584400e-9 * ta * vel3 * d_tr * pa_pr + \
        2.20609296e-7 * vel4 * d_tr * pa_pr + \
        -7.32469180e-4 * d_tr2 * pa_pr + \
        -1.87381964e-5 * ta * d_tr2 * pa_pr + \
        4.80925239e-6 * ta2 * d_tr2 * pa_pr + \
        -8.75492040e-8 * ta3 * d_tr2 * pa_pr + \
        2.77862930e-5 * vel * d_tr2 * pa_pr + \
        -5.06004592e-6 * ta * vel * d_tr2 * pa_pr + \
        1.14325367e-7 * ta2 * vel * d_tr2 * pa_pr + \
        2.53016723e-6 * vel2 * d_tr2 * pa_pr + \
        -1.72857035e-8 * ta * vel2 * d_tr2 * pa_pr + \
        -3.95079398e-8 * vel3 * d_tr2 * pa_pr + \
        -3.59413173e-7 * d_tr3 * pa_pr + \
        7.04388046e-7 * ta * d_tr3 * pa_pr + \
        -1.89309167e-8 * ta2 * d_tr3 * pa_pr + \
        -4.79768731e-7 * vel * d_tr3 * pa_pr + \
        7.96079978e-9 * ta * vel * d_tr3 * pa_pr + \
        1.62897058e-9 * vel2 * d_tr3 * pa_pr + \
        3.94367674e-8 * d_tr4 * pa_pr + \
        -1.18566247e-9 * ta * d_tr4 * pa_pr + \
        3.34678041e-10 * vel * d_tr4 * pa_pr + \
        -1.15606447e-10 * d_tr5 * pa_pr + \
        -2.80626406 * pa_pr2 + \
        0.548712484 * ta * pa_pr2 + \
        -0.00399428410 * ta2 * pa_pr2 + \
        -9.54009191e-4 * ta3 * pa_pr2 + \
        1.93090978e-5 * ta4 * pa_pr2 + \
        -0.308806365 * vel * pa_pr2 + \
        0.0116952364 * ta * vel * pa_pr2 + \
        4.95271903e-4 * ta2 * vel * pa_pr2 + \
        -1.90710882e-5 * ta3 * vel * pa_pr2 + \
        0.00210787756 * vel2 * pa_pr2 + \
        -6.98445738e-4 * ta * vel2 * pa_pr2 + \
        2.30109073e-5 * ta2 * vel2 * pa_pr2 + \
        4.17856590e-4 * vel3 * pa_pr2 + \
        -1.27043871e-5 * ta * vel3 * pa_pr2 + \
        -3.04620472e-6 * vel4 * pa_pr2 + \
        0.0514507424 * d_tr * pa_pr2 + \
        -0.00432510997 * ta * d_tr * pa_pr2 + \
        8.99281156e-5 * ta2 * d_tr * pa_pr2 + \
        -7.14663943e-7 * ta3 * d_tr * pa_pr2 + \
        -2.66016305e-4 * vel * d_tr * pa_pr2 + \
        2.63789586e-4 * ta * vel * d_tr * pa_pr2 + \
        -7.01199003e-6 * ta2 * vel * d_tr * pa_pr2 + \
        -1.06823306e-4 * vel2 * d_tr * pa_pr2 + \
        3.61341136e-6 * ta * vel2 * d_tr * pa_pr2 + \
        2.29748967e-7 * vel3 * d_tr * pa_pr2 + \
        3.04788893e-4 * d_tr2 * pa_pr2 + \
        -6.42070836e-5 * ta * d_tr2 * pa_pr2 + \
        1.16257971e-6 * ta2 * d_tr2 * pa_pr2 + \
        7.68023384e-6 * vel * d_tr2 * pa_pr2 + \
        -5.47446896e-7 * ta * vel * d_tr2 * pa_pr2 + \
        -3.59937910e-8 * vel2 * d_tr2 * pa_pr2 + \
        -4.36497725e-6 * d_tr3 * pa_pr2 + \
        1.68737969e-7 * ta * d_tr3 * pa_pr2 + \
        2.67489271e-8 * vel * d_tr3 * pa_pr2 + \
        3.23926897e-9 * d_tr4 * pa_pr2 + \
        -0.0353874123 * pa_pr3 + \
        -0.221201190 * ta * pa_pr3 + \
        0.0155126038 * ta2 * pa_pr3 + \
        -2.63917279e-4 * ta3 * pa_pr3 + \
        0.0453433455 * vel * pa_pr3 + \
        -0.00432943862 * ta * vel * pa_pr3 + \
        1.45389826e-4 * ta2 * vel * pa_pr3 + \
        2.17508610e-4 * vel2 * pa_pr3 + \
        -6.66724702e-5 * ta * vel2 * pa_pr3 + \
        3.33217140e-5 * vel3 * pa_pr3 + \
        -0.00226921615 * d_tr * pa_pr3 + \
        3.80261982e-4 * ta * d_tr * pa_pr3 + \
        -5.45314314e-9 * ta2 * d_tr * pa_pr3 + \
        -7.96355448e-4 * vel * d_tr * pa_pr3 + \
        2.53458034e-5 * ta * vel * d_tr * pa_pr3 + \
        -6.31223658e-6 * vel2 * d_tr * pa_pr3 + \
        3.02122035e-4 * d_tr2 * pa_pr3 + \
        -4.77403547e-6 * ta * d_tr2 * pa_pr3 + \
        1.73825715e-6 * vel * d_tr2 * pa_pr3 + \
        -4.09087898e-7 * d_tr3 * pa_pr3 + \
        0.614155345 * pa_pr4 + \
        -0.0616755931 * ta * pa_pr4 + \
        0.00133374846 * ta2 * pa_pr4 + \
        0.00355375387 * vel * pa_pr4 + \
        -5.13027851e-4 * ta * vel * pa_pr4 + \
        1.02449757e-4 * vel2 * pa_pr4 + \
        -0.00148526421 * d_tr * pa_pr4 + \
        -4.11469183e-5 * ta * d_tr * pa_pr4 + \
        -6.80434415e-6 * vel * d_tr * pa_pr4 + \
        -9.77675906e-6 * d_tr2 * pa_pr4 + \
        0.0882773108 * pa_pr5 + \
        -0.00301859306 * ta * pa_pr5 + \
        0.00104452989 * vel * pa_pr5 + \
        2.47090539e-4 * d_tr * pa_pr5 + \
        0.00148348065 * pa_pr6
    return utci_approx



def hottest_day_date(_air_temp):
    hottest_day_value = max(_air_temp)
    hottest_day_index = _air_temp.index(hottest_day_value) +1

    day_of_year = hottest_day_index // 24 + 1
    hour_of_day = hottest_day_index % 24

    year = datetime.datetime.now().year
    date = datetime.datetime(year, 1, 1) + datetime.timedelta(
        days=day_of_year - 1,
        hours=hour_of_day
    )
    month = date.month
    day = date.day
    formatted_date = date.strftime("%d %b %H:%M")
    return formatted_date, month, day, hottest_day_index 
    
def average_value(list):
    return gh.Average(list)
    

if run:
    _air_temp, _wind_vel_, _rel_humid = epw_open(EPW)
    date_in_hottestDay, month, day, hottestDay_index = hottest_day_date(_air_temp)
    Epw_year = "Data provison from {} to {}" .format(max(model_year), min(model_year))
    avg_temp = max(_air_temp)
    avg_humidity = gh.Average(relative_humidity) 
    wind_velocity_on_hottestDay = wind_speed[hottestDay_index]   
    
    utci_list = [
        universal_thermal_climate_index(ta_i, v_i, rh_i)
        for ta_i, v_i, rh_i in zip(_air_temp, _wind_vel_, _rel_humid)
    ]
    hottest_utci = max(utci_list)
    period = "{}/{} to {}/{} between 0 and 23 @1" .format(month, day, month, day)