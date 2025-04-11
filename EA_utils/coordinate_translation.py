def utm_to_latlon(eastings, northings, zone_number, hemisphere):
    #usage example:
    # import EA_utils.coordinate_translation as ct
    #####################################################
    # eastings = [852934.6078730285]
    # northings = [6097546.724721148]
    # zone_number = 32
    # hemisphere = "N"
    #####################################################
    # _ = ct.utm_to_latlon(eastings, northings, zone_number, hemisphere)
    import utm

    latitudes = [[]]*len(eastings)
    longitudes = [[]]*len(eastings)
    for location_idx, easting in enumerate(eastings):
        northing = northings[location_idx]
        if hemisphere == "N":
            (latitudes[location_idx],longitudes[location_idx]) = utm.to_latlon(easting, northing, zone_number, northern=True)
        elif hemisphere == "S":
            (latitudes[location_idx],longitudes[location_idx]) = utm.to_latlon(easting, northing, zone_number, northern=False)
        else:
            raise ValueError("hemisphere must be ""N"" or ""S""")

        print(f'location #{location_idx+1}:')
        print(f"latitude = {latitudes[location_idx]}")
        print(f"longitude = {longitudes[location_idx]}")
        if location_idx<len(eastings)-1: print(f"\n")

    return latitudes, longitudes

def latlon_to_utm(latitudes, longitudes, zone_number=None):
    #usage example:
    # import EA_utils.coordinate_translation as ct
    #####################################################
    # latitudes = [54.494569]
    # longitudes = [5.775285]
    # zone_number = 32 #define or leave empty to derive and output
    #####################################################
    # ct.latlon_to_utm(latitudes, longitudes, zone_number=None)
    from pyproj import Proj
    import utm
    
    eastings = [[]]*len(latitudes)
    northings = [[]]*len(latitudes)
    for location_idx, latitude in enumerate(latitudes):
        longitude = longitudes[location_idx]
        if zone_number:
            myProj = Proj(proj='utm', zone=f"{zone_number}", ellps='WGS84')
            eastings[location_idx], northings[location_idx] = myProj(longitude, latitude)
            zone_letter = []
        else:
            (eastings[location_idx], northings[location_idx], zone_number, zone_letter) = utm.from_latlon(latitude, longitude)
        print(f'location #{location_idx+1}:')
        print(f"easting = {eastings[location_idx]}")
        print(f"northing = {northings[location_idx]}")
        print(f"zone number = {zone_number}")
        if zone_letter: print(f"zone letter = {zone_letter}")
        if location_idx<len(latitudes)-1: print(f"\n")
    
    return eastings, northings, zone_number, zone_letter