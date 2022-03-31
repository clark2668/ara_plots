from stringprep import in_table_c21_c22
import numpy as np
FTtoM = 0.3048
MtoFT = 1./FTtoM

a2_loc = np.array( [-0.5986437007,  0.8010086064, -0.000970503])


def get_angle(thevector):
    return np.rad2deg(np.arctan2(thevector[1], thevector[0]))

def transform_map_to_global(input_azi):
    angle_from_zero = np.rad2deg(np.arctan2(a2_loc[1], a2_loc[0]))
    output_azi = input_azi + angle_from_zero
    if output_azi > 180.:
        output_azi -= 360
    return output_azi

# double TransformMapPeakToGlobalFrame(double inputAzimuth){

#     double TBx[3] = {-0.5986468972,0.8010128835,-0.000399302}; //this is from the ARA coordinates document, and contains the easting, northing, and upping components of the TB origin relative to the array (and therefore global) origin                                                                                                                                                                                                                                             
#     double angle_from_0 = TMath::ATan2(TBx[1],TBx[0]) * TMath::RadToDeg(); //I can compute this as the angle from zero because ATan2 is intelligent, and knows I mean the angle from zero!                                                
#     double outputAzimuth = inputAzimuth + angle_from_0; //correct this                                                                                                                                                                    
#     if(outputAzimuth > 180.) outputAzimuth -=360.; //this accounts for domain correcting the answer (ie, spinning from -180 -> 180, to 0->360                                                                                             
#     return outputAzimuth; //return the answer                                                                                                                                                                                             

# }

# double TransformGlobalFrameToMapFrame(double inputAzimuth){

#     double TBx[3] = {-0.5986468972,0.8010128835,-0.000399302}; //this is from the ARA coordinates document, and contains the easting, northing, and upping components of the TB origin relative to the array (and therefore global) origin

#     double angle_from_0 = TMath::ATan2(TBx[1],TBx[0]) * TMath::RadToDeg(); 
#     double outputAzimuth = inputAzimuth - angle_from_0;
#     if(outputAzimuth < 180.) outputAzimuth+=360.;
#     return outputAzimuth;

# }

def get_conversion(mode):
    conversion = 1.

    if mode=='meters':
        conversion = FTtoM
    elif mode=='km':
        conversion = FTtoM/1000.
    elif mode=='ft':
        conversion = 1
    return conversion

def get_icecube_string_locations(mode='meters'):

    conversion = get_conversion(mode)

    data = np.genfromtxt("icecube_dom60.csv", delimiter=',',skip_header=1, names=['x','y','z','string','dom'])

    xs = data['x'] * MtoFT * conversion
    ys = data['y'] * MtoFT * conversion
    zs = data['z'] * MtoFT * conversion
    strings = data['string']
    doms = data['dom']
    return xs, ys


def get_thing(thing, mode='meters'):

    conversion = get_conversion(mode)

    if thing == 'A2':
        return np.array([35528.80, 45382.35]) * conversion
    elif thing == 'I3Center':
        return np.array([46500.,52200.]) * conversion
    elif thing =='SPIce':
        return np.array([42560.,48790.]) * conversion
    elif thing =='IC1':
        # we enter all values as feet by default
        ic22 = np.array([-256.14,-521.08]) * MtoFT * conversion
        i3center = get_thing('I3Center', mode)
        ic22[0]+=i3center[0]
        ic22[1]+=i3center[1]
        return ic22
    elif thing =='IC22':
        # we enter all values as feet by default
        ic22 = np.array([-492.43,-230.16]) * MtoFT * conversion
        i3center = get_thing('I3Center', mode)
        ic22[0]+=i3center[0]
        ic22[1]+=i3center[1]
        return ic22
    elif thing =='RTP':
        return np.array([24000.78+22400,-1728.47 + 53900]) * conversion
    # elif thing == 'IC1':
    #     return 