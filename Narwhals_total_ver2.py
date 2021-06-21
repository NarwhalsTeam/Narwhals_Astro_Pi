# -*- coding: utf-8 -*-


import math
import numpy as np
from ephem import readtle, degree
from datetime import datetime
import csv

import glob
import matplotlib.pyplot as plt


import os, unittest
from datetime import date





# Geomag WMM model
# the code was taken from:
# https://www.ngdc.noaa.gov/geomag/WMM/soft.shtml
# https://www.ngdc.noaa.gov/geomag/WMM/thirdpartycontributions.shtml
# https://pypi.org/project/geomag/


class GeoMag:

    def GeoMag(self, dlat, dlon, h=0, time=date.today()): # latitude (decimal degrees), longitude (decimal degrees), altitude (feet), date
        #time = date('Y') + date('z')/365
        time = time.year+((time - date(time.year,1,1)).days/365.0)
        
        # alt = h/3280.8399

        alt = h
        otime = oalt = olat = olon = -1000.0

        dt = time - self.epoch
        glat = dlat
        glon = dlon
        rlat = math.radians(glat)
        rlon = math.radians(glon)
        srlon = math.sin(rlon)
        srlat = math.sin(rlat)
        crlon = math.cos(rlon)
        crlat = math.cos(rlat)
        srlat2 = srlat*srlat
        crlat2 = crlat*crlat
        self.sp[1] = srlon
        self.cp[1] = crlon

        #/* CONVERT FROM GEODETIC COORDS. TO SPHERICAL COORDS. */
        if (alt != oalt or glat != olat):
            q = math.sqrt(self.a2-self.c2*srlat2)
            q1 = alt*q
            q2 = ((q1+self.a2)/(q1+self.b2))*((q1+self.a2)/(q1+self.b2))
            ct = srlat/math.sqrt(q2*crlat2+srlat2)
            st = math.sqrt(1.0-(ct*ct))
            r2 = (alt*alt)+2.0*q1+(self.a4-self.c4*srlat2)/(q*q)
            r = math.sqrt(r2)
            d = math.sqrt(self.a2*crlat2+self.b2*srlat2)
            ca = (alt+d)/r
            sa = self.c2*crlat*srlat/(r*d)

        if (glon != olon):
            for m in range(2,self.maxord+1):
                self.sp[m] = self.sp[1]*self.cp[m-1]+self.cp[1]*self.sp[m-1]
                self.cp[m] = self.cp[1]*self.cp[m-1]-self.sp[1]*self.sp[m-1]

        aor = self.re/r
        ar = aor*aor
        br = bt = bp = bpp = 0.0
        for n in range(1,self.maxord+1):
            ar = ar*aor
            
            #for (m=0,D3=1,D4=(n+m+D3)/D3;D4>0;D4--,m+=D3):
            m=0
            D3=1
            #D4=(n+m+D3)/D3
            D4=(n+m+1)
            while D4>0:

        # /*
                # COMPUTE UNNORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS
                # AND DERIVATIVES VIA RECURSION RELATIONS
        # */
                if (alt != oalt or glat != olat):
                    if (n == m):
                        self.p[m][n] = st * self.p[m-1][n-1]
                        self.dp[m][n] = st*self.dp[m-1][n-1]+ct*self.p[m-1][n-1]

                    elif (n == 1 and m == 0):
                        self.p[m][n] = ct*self.p[m][n-1]
                        self.dp[m][n] = ct*self.dp[m][n-1]-st*self.p[m][n-1]

                    elif (n > 1 and n != m):
                        if (m > n-2):
                            self.p[m][n-2] = 0
                        if (m > n-2):
                            self.dp[m][n-2] = 0.0
                        self.p[m][n] = ct*self.p[m][n-1]-self.k[m][n]*self.p[m][n-2]
                        self.dp[m][n] = ct*self.dp[m][n-1] - st*self.p[m][n-1]-self.k[m][n]*self.dp[m][n-2]

        # /*
                # TIME ADJUST THE GAUSS COEFFICIENTS
        # */
                if (time != otime):
                    self.tc[m][n] = self.c[m][n]+dt*self.cd[m][n]
                    if (m != 0):
                        self.tc[n][m-1] = self.c[n][m-1]+dt*self.cd[n][m-1]

        # /*
                # ACCUMULATE TERMS OF THE SPHERICAL HARMONIC EXPANSIONS
        # */
                par = ar*self.p[m][n]
                
                if (m == 0):
                    temp1 = self.tc[m][n]*self.cp[m]
                    temp2 = self.tc[m][n]*self.sp[m]
                else:
                    temp1 = self.tc[m][n]*self.cp[m]+self.tc[n][m-1]*self.sp[m]
                    temp2 = self.tc[m][n]*self.sp[m]-self.tc[n][m-1]*self.cp[m]

                bt = bt-ar*temp1*self.dp[m][n]
                bp = bp + (self.fm[m] * temp2 * par)
                br = br + (self.fn[n] * temp1 * par)
        # /*
                    # SPECIAL CASE:  NORTH/SOUTH GEOGRAPHIC POLES
        # */
                if (st == 0.0 and m == 1):
                    if (n == 1):
                        self.pp[n] = self.pp[n-1]
                    else:
                        self.pp[n] = ct*self.pp[n-1]-self.k[m][n]*self.pp[n-2]
                    parp = ar*self.pp[n]
                    bpp = bpp + (self.fm[m]*temp2*parp)
                    
                D4=D4-1
                m=m+1

        if (st == 0.0):
            bp = bpp
        else:
            bp = bp/st
        # /*
            # ROTATE MAGNETIC VECTOR COMPONENTS FROM SPHERICAL TO
            # GEODETIC COORDINATES
        # */
        bx = -bt*ca-br*sa
        by = bp
        bz = bt*sa-br*ca
        # /*
            # COMPUTE DECLINATION (DEC), INCLINATION (DIP) AND
            # TOTAL INTENSITY (TI)
        # */
        bh = math.sqrt((bx*bx)+(by*by))
        ti = math.sqrt((bh*bh)+(bz*bz))
        dec = math.degrees(math.atan2(by,bx))
        dip = math.degrees(math.atan2(bz,bh))
        # /*
            # COMPUTE MAGNETIC GRID VARIATION IF THE CURRENT
            # GEODETIC POSITION IS IN THE ARCTIC OR ANTARCTIC
            # (I.E. GLAT > +55 DEGREES OR GLAT < -55 DEGREES)

            # OTHERWISE, SET MAGNETIC GRID VARIATION TO -999.0
        # */
        gv = -999.0
        if (math.fabs(glat) >= 55.):
            if (glat > 0.0 and glon >= 0.0):
                gv = dec-glon
            if (glat > 0.0 and glon < 0.0):
                gv = dec+math.fabs(glon);
            if (glat < 0.0 and glon >= 0.0):
                gv = dec+glon
            if (glat < 0.0 and glon < 0.0):
                gv = dec-math.fabs(glon)
            if (gv > +180.0):
                gv = gv - 360.0
            if (gv < -180.0):
                gv = gv + 360.0

        otime = time
        oalt = alt
        olat = glat
        olon = glon

        class RetObj:
            pass
        retobj = RetObj()
        retobj.dec = dec
        retobj.dip = dip
        retobj.ti = ti
        retobj.bh = bh
        retobj.bx = bx
        retobj.by = by
        retobj.bz = bz
        retobj.lat = dlat
        retobj.lon = dlon
        retobj.alt = h
        retobj.time = time

        return retobj

    def __init__(self, wmm_filename=None):
        if not wmm_filename:
            wmm_filename = os.path.join(os.path.dirname(__file__), 'WMM.COF')
        wmm=[]
        with open(wmm_filename) as wmm_file:
            for line in wmm_file:
                linevals = line.strip().split()
                if len(linevals) == 3:
                    self.epoch = float(linevals[0])
                    self.model = linevals[1]
                    self.modeldate = linevals[2]
                elif len(linevals) == 6:
                    linedict = {'n': int(float(linevals[0])),
                    'm': int(float(linevals[1])),
                    'gnm': float(linevals[2]),
                    'hnm': float(linevals[3]),
                    'dgnm': float(linevals[4]),
                    'dhnm': float(linevals[5])}
                    wmm.append(linedict)

        z = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
        self.maxord = self.maxdeg = 12
        self.tc = [z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13]]
        self.sp = z[0:14]
        self.cp = z[0:14]
        self.cp[0] = 1.0
        self.pp = z[0:13]
        self.pp[0] = 1.0
        self.p = [z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14]]
        self.p[0][0] = 1.0
        self.dp = [z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13]]
        self.a = 6378.137
        self.b = 6356.7523142
        self.re = 6371.2
        self.a2 = self.a*self.a
        self.b2 = self.b*self.b
        self.c2 = self.a2-self.b2
        self.a4 = self.a2*self.a2
        self.b4 = self.b2*self.b2
        self.c4 = self.a4 - self.b4

        self.c = [z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14]]
        self.cd = [z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14],z[0:14]]
        
        for wmmnm in wmm:
            m = wmmnm['m']
            n = wmmnm['n']
            gnm = wmmnm['gnm']
            hnm = wmmnm['hnm']
            dgnm = wmmnm['dgnm']
            dhnm = wmmnm['dhnm']
            if (m <= n):
                self.c[m][n] = gnm
                self.cd[m][n] = dgnm
                if (m != 0):
                    self.c[n][m-1] = hnm
                    self.cd[n][m-1] = dhnm

        #/* CONVERT SCHMIDT NORMALIZED GAUSS COEFFICIENTS TO UNNORMALIZED */
        self.snorm = [z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13]]
        self.snorm[0][0] = 1.0
        self.k = [z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13],z[0:13]]
        self.k[1][1] = 0.0
        self.fn = [0.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0]
        self.fm = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0]
        for n in range(1,self.maxord+1):
            self.snorm[0][n] = self.snorm[0][n-1]*(2.0*n-1)/n
            j=2.0
            #for (m=0,D1=1,D2=(n-m+D1)/D1;D2>0;D2--,m+=D1):
            m=0
            D1=1
            D2=(n-m+D1)/D1
            while (D2 > 0):
                self.k[m][n] = (((n-1)*(n-1))-(m*m))/((2.0*n-1)*(2.0*n-3.0))
                if (m > 0):
                    flnmj = ((n-m+1.0)*j)/(n+m)
                    self.snorm[m][n] = self.snorm[m-1][n]*math.sqrt(flnmj)
                    j = 1.0
                    self.c[n][m-1] = self.snorm[m][n]*self.c[n][m-1]
                    self.cd[n][m-1] = self.snorm[m][n]*self.cd[n][m-1]
                self.c[m][n] = self.snorm[m][n]*self.c[m][n]
                self.cd[m][n] = self.snorm[m][n]*self.cd[m][n]
                D2=D2-1
                m=m+D1

class GeoMagTest(unittest.TestCase):

    d1=date(2015,1,1)
    d2=date(2017,7,2)
    
    test_values = (
        # date, alt, lat, lon, var
        (d1, 0, 80, 0,  -3.85),
        (d1, 0, 0, 120, 0.57),
        (d1, 0, -80, 240,  69.81),
        (d1, 328083.99, 80, 0, -4.27),
        (d1, 328083.99, 0, 120, 0.56),
        (d1, 328083.99, -80, 240, 69.22),
        (d2, 0, 80, 0, -2.75),
        (d2, 0, 0, 120, 0.32),
        (d2, 0, -80, 240, 69.58),
        (d2, 328083.99, 80, 0, -3.17),
        (d2, 328083.99, 0, 120, 0.32),
        (d2, 328083.99, -80, 240, 69.00),
    )
    
    def test_declination(self):
        gm = GeoMag()
        for values in self.test_values:
            calcval=gm.GeoMag(values[2], values[3], values[1], values[0])
            self.assertAlmostEqual(values[4], calcval.dec, 2, 'Expected %s, result %s' % (values[4], calcval.dec))





""" Latest TLE data for ISS location"""
name = "ISS (ZARYA)"  


#ED 22.04.2021

line1 = "1 25544U 98067A   21112.22487556  .00001865  00000-0  42113-4 0  9991"
line2 = "2 25544  51.6450 259.5242 0002599 262.0008 158.2898 15.48917530279928"        

# IZZY
#line1 = "1 25544U 98067A   16069.57438447  .00010945  00000-0  17329-3 0  9998"
#line2 = "2 25544  51.6422 194.8005 0001583 253.9754 212.9290 15.53974450989487"           


""" Coordinates of the South Magnetic Pole wich is close to the geographic North Pole """
north_lat = 80.65
north_long = -72.68

def vincenty(coordA, coordB):
    """
    Calculates the angle between 2 points on Earth with coordinates coordA and coordB
    respectively. Typically, in our project, the second pair represents the coordinates 
    of the South Magnetic Pole
    """
    """Constants"""
    RE = 6378137.0                                 # radius at equator in meters (WGS-84)
    FLAT = 1 / 298.257223563                       # flattening of the ellipsoid (WGS-84)
    B = (1 - FLAT) * RE                

    latitudeA, longitudeA = coordA                 
    latitudeB, longitudeB = coordB
    u_1 = math.atan((1 - FLAT) * math.tan(math.radians(latitudeA)))
    u_2 = math.atan((1 - FLAT) * math.tan(math.radians(latitudeB)))

    L = math.radians(longitudeB - longitudeA)

    Lambda = L                                # set initial value of lambda to L
    sin_u1 = math.sin(u_1)
    cos_u1 = math.cos(u_1)
    sin_u2 = math.sin(u_2)
    cos_u2 = math.cos(u_2)
    
    """ Short-circuit coincident point """
    if latitudeA == latitudeB and longitudeA == longitudeB:
        return 0.0
    
    """ Convergence iterations """
    for _ in range (0, 300):
       cos_lambda = math.cos(Lambda)
       sin_lambda = math.sin(Lambda)
       sin_sigma = math.sqrt((cos_u2 * math.sin(Lambda)) ** 2 + 
                    (cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lambda) ** 2)
       if sin_sigma == 0:
           return 0.0  
       cos_sigma = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_lambda
       sigma = math.atan2(sin_sigma, cos_sigma)
          
       sin_alpha = (cos_u1 * cos_u2 * sin_lambda) / sin_sigma
       cos_sq_alpha = 1 - sin_alpha ** 2
       try:
           cos2_sigma_m = cos_sigma - ((2 * sin_u1 * sin_u2) / cos_sq_alpha)
       except ZeroDivisionError:
           cos2_sigma_m = 0
            
       C = (FLAT / 16) * cos_sq_alpha * (4 + FLAT *(4 - 3 * cos_sq_alpha))
       Lambda_prev = Lambda
       Lambda = L + (1 - C) * FLAT * sin_alpha * (sigma + C * sin_sigma * 
                (cos2_sigma_m + C * cos_sigma * (-1 + 2 * cos2_sigma_m ** 2)))
        
       """ Successful convergence """
       diff = abs(Lambda_prev - Lambda)
       if diff <= 10**-12:
           break
        
    u_sq = cos_sq_alpha * ((RE ** 2 - B ** 2) / B ** 2)
    X = 1 + (u_sq / 16384) * (4096 + u_sq * (-768 + u_sq * (320 - 175 * u_sq)))
    Y = (u_sq / 1024) * (256 + u_sq * (-128 + u_sq * (74 - 47 * u_sq)))
    delta_sig = Y * sin_sigma * (cos2_sigma_m + 0.25 * Y * (cos_sigma * 
                (-1 + 2 * cos2_sigma_m ** 2) - (1 / 6) * Y * cos2_sigma_m *
                (-3 + 4 * sin_sigma ** 2) * (-3 + 4 * cos2_sigma_m ** 2)))

    distance = B * X * (sigma - delta_sig) # Distance in meters from point A to point B 

    return distance / RE  # The angle Theta from point A to point B in radians 

""" Constants """
B0 = 3.12 * 1e-5
RE = 6370
r = RE + 408

    
def calculate_magnetic_field_t(theta):
    """
    Calculate the ideal magnetic field for a certain position of the ISS, knowing the angle
    *theta*, the polar angle,  between the ISS location and the South Magnetic Pole.
    Formula:
    |B| = B0 * (RE / r) ^ 3 * sqrt(1 + 3 * cos^2(theta))
    
    Where:
    B0 is a constant = 3.12 * 10^-5
    RE is the Radius of the Earth
    r is the distance from the ISS to the center of the Earth
    """
    

    """ Formula """
    B = B0 * (RE / r) ** 3
    B = B * math.sqrt(1 + 3 * (math.cos(theta)) ** 2)
    B = B * 10**6     # Transform B from tesla to microtesla
    return B

def calculate_magnetic_field_z(theta):
    B = B0 * (RE / r) ** 3
    B = B *2 * math.cos(theta)
    B = B * 10**6 
    return B

def calculate_magnetic_field_x(theta):
    B = B0 * (RE / r) ** 3
    B = B * math.sin(theta)
    B = B * 10**6 
    return B    
    

# MAIN PROGRAM


gm = GeoMag ("WMM.COF")

with open('total_interval_sec_prelucrat.csv', 'w',newline='') as f_output:
    csv_writer_f_output = csv.writer(f_output)

    with open('total_interval_sec.csv', 'r',newline='') as f_input:
        csv_reader_f_input = csv.reader(f_input)
        line_count = 0
        for row in csv_reader_f_input:
            if line_count == 0:
                #print(f' Header  {", ".join(row)}')
                row_work=['row','time_REC', 'pitch_rec', 'roll_rec', 'yaw_rec','beta','lat_rec_r', 'long_rec_r','lat_r', 'long_r','Theta_rad', 'B_x', 'B_y', 'B_z', 'B_t' ,'Bth_x', 'Bth_y' , 'Bth_z', 'Bth_t', 'Bwmm_x', 'Bwmm_y' , 'Bwmm_z', 'Bwmm_t','Bwmm1_x', 'Bwmm1_y' , 'Bwmm1_z', 'Bwmm2_x', 'Bwmm2_y' , 'Bwmm2_z','Bwmm3_x', 'Bwmm3_y' , 'Bwmm3_z','Bwmm3_t', 'pitch_api', 'roll_api', 'yaw_api','Biss_x', 'Biss_API_y' , 'Biss_API_z', 'Biss_API_t','B_p_x', 'B_p_y', 'B_p_z', 'B_p_t', 'pitch_iss','roll_iss', 'yaw_iss','Bwmm_p_x','Bwmm_p_y' , 'Bwmm_p_z', 'Bwmm_p_t', 'Biss_VEC_x', 'Biss_VEC_y' , 'Biss_VEC_z', 'Biss_VEC_t', 'alpha', 'vt_Bwmm_paralel', 'vs_Bwmm_paralel', 'I_indus', 'mom_mag_indus', 'Forta_wmm_indus', 'mom_F_wmm_indus', 'mom_mag_sursa', 'Forta_wmm_sursa', 'mom_F_wmm_sursa' ]
                # header= row+row_work
                header = row_work
                csv_writer_f_output.writerow(header)
                line_count += 1
                beta=0
            else:
               
                # Latitude si longitude using the time and TLE data 
                datetime_calc=row[0]
                now_time = datetime.now()
                                             
              
                iss = readtle(name, line1, line2)
                iss.compute(datetime_calc)
                lat = iss.sublat/degree
                long = iss.sublong/degree
                
                lat_rec_r = float(row[7])  * 3.14/180
                long_rec_r = float(row[8]) * 3.14/180
                
                               
                # Recorded EULER angles
                pitch_rec = float(row[15])
                roll_rec = float(row[16])
                yaw_rec = float(row[17])
                
            
                 
                pitch_rec = math.radians(pitch_rec)
                roll_rec = math.radians(roll_rec) 
                yaw_rec =  math.radians(yaw_rec)
                
                # computing of beta, the angle between ISS’ velocity and the geographic North direction
                lat_r = lat * 3.14/180
                long_r = long * 3.14/180
                
                if line_count == 1: 
                    lat_r_ant = lat_r
                    long_r_ant = long_r
                
                d_lat_r = lat_r - lat_r_ant
                d_long_r = long_r - long_r_ant
                
                lat_r_ant = lat_r
                long_r_ant = long_r
                
                    
                if d_long_r != 0 and d_long_r < 0.2 and d_long_r > -0.2:
                    raport =  d_lat_r/d_long_r
                    if raport > 0 :
                        beta = math.atan(1/raport)
                    else :
                        beta = 3.14 + math.atan(1/ raport)
                else:
                    beta = beta  
                
               
                                                          
                # Measured magnetic field in SR API
                B_x = float(row[10])
                B_y = float(row[11])
                B_z = float(row[12])
                
                B_t = math.sqrt(B_x**2 + B_y**2  + B_z**2)
                
                magnetic_field_Matrix = np.matrix([
                [B_x],
                [B_y],
                [B_z]
                ])
                
               
                # Magetic field : theorethical dipol  model
                Theta = vincenty([lat, long], [north_lat, north_long])
                Theta_rad = (Theta * 3.14) / 180
                Bth_x = calculate_magnetic_field_x(Theta)
                Bth_z = calculate_magnetic_field_z(Theta)
                Bth_y = 0
                Bth_t = calculate_magnetic_field_t(Theta)
              
                # WMM
                date_time_obj =  datetime.strptime(datetime_calc, '%Y-%m-%d %H:%M:%S.%f') 
                work_time = date_time_obj.date()
               
                
                               
                            
                # Magetic field : WMM
                mag = gm.GeoMag(lat,long,408, work_time)
                Bwmm_x = mag.bx / 1000
                Bwmm_y = mag.by / 1000
                Bwmm_z = mag.bz / 1000  
                Bwmm_t = math.sqrt(Bwmm_x**2 + Bwmm_y**2  + Bwmm_z**2)   
                
                                  
                magnetic_field_Matrix_wmm = np.matrix([
                [Bwmm_x],
                [Bwmm_y],
                [Bwmm_z]
                ])
                
                      
                           
                
               
                # SR API -  Transformation WMM NEC --> API
                
                # Measured angles
                pitch_api = pitch_rec
                roll_api = roll_rec
                yaw_api =  yaw_rec
                
                # Direct Euler matrix
                yawMatrix_NEC_API= np.matrix([
                [math.cos(yaw_api), math.sin(yaw_api), 0],
                [-math.sin(yaw_api), math.cos(yaw_api), 0],
                [0, 0, 1]
                ])

                pitchMatrix_NEC_API= np.matrix([
                [math.cos(pitch_api), 0, -math.sin(pitch_api)],
                [0, 1, 0],
                [math.sin(pitch_api), 0, math.cos(pitch_api)]
                ])

                rollMatrix_NEC_API = np.matrix([
                [1, 0, 0],
                [0, math.cos(roll_api), math.sin(roll_api)],
                [0, -math.sin(roll_api), math.cos(roll_api)]
                ])
                
                R_NEC_API = rollMatrix_NEC_API * pitchMatrix_NEC_API *yawMatrix_NEC_API
                 
             
                # Magnetic field : processing
                magnetic_field_p_Matrix_wmm = R_NEC_API *  magnetic_field_Matrix_wmm
                
                Bwmm1_x = magnetic_field_p_Matrix_wmm[0,0]
                Bwmm1_y = magnetic_field_p_Matrix_wmm[1,0]
                Bwmm1_z = magnetic_field_p_Matrix_wmm[2,0]
                
                
                # Filter for pitch and roll, YAW = Beta 
                pitch_rec_flt = 0.1
                roll_rec_flt = 1.8
                
               
                pitch_api = pitch_rec_flt 
                roll_api = roll_rec_flt
                yaw_api =  beta
                
                # Direct Euler matrix
                yawMatrix_NEC_API= np.matrix([
                [math.cos(yaw_api), math.sin(yaw_api), 0],
                [-math.sin(yaw_api), math.cos(yaw_api), 0],
                [0, 0, 1]
                ])

                pitchMatrix_NEC_API= np.matrix([
                [math.cos(pitch_api), 0, -math.sin(pitch_api)],
                [0, 1, 0],
                [math.sin(pitch_api), 0, math.cos(pitch_api)]
                ])

                rollMatrix_NEC_API = np.matrix([
                [1, 0, 0],
                [0, math.cos(roll_api), math.sin(roll_api)],
                [0, -math.sin(roll_api), math.cos(roll_api)]
                ])
                
                R_NEC_API = rollMatrix_NEC_API * pitchMatrix_NEC_API *yawMatrix_NEC_API
                 
             
                # Magnetic field : processing
                magnetic_field_p_Matrix_wmm = R_NEC_API *  magnetic_field_Matrix_wmm
                
                Bwmm2_x = magnetic_field_p_Matrix_wmm[0,0]
                Bwmm2_y = magnetic_field_p_Matrix_wmm[1,0]
                Bwmm2_z = magnetic_field_p_Matrix_wmm[2,0]
                
               
                
                
                #  Rotatin to match WMM and the measured field
                pitch_api =  0.9 + pitch_rec_flt
                roll_api =   3.14 + roll_rec_flt
                yaw_api =  beta  + 1.57  
                
                # Direct Euler matrix
                yawMatrix_NEC_API= np.matrix([
                [math.cos(yaw_api), math.sin(yaw_api), 0],
                [-math.sin(yaw_api), math.cos(yaw_api), 0],
                [0, 0, 1]
                ])

                pitchMatrix_NEC_API= np.matrix([
                [math.cos(pitch_api), 0, -math.sin(pitch_api)],
                [0, 1, 0],
                [math.sin(pitch_api), 0, math.cos(pitch_api)]
                ])

                rollMatrix_NEC_API = np.matrix([
                [1, 0, 0],
                [0, math.cos(roll_api), math.sin(roll_api)],
                [0, -math.sin(roll_api), math.cos(roll_api)]
                ])
                
                R_NEC_API = rollMatrix_NEC_API * pitchMatrix_NEC_API *yawMatrix_NEC_API
                 
  
            
                # Magnetic field : processing
                magnetic_field_p_Matrix_wmm = R_NEC_API *  magnetic_field_Matrix_wmm
                
                Bwmm3_x = magnetic_field_p_Matrix_wmm[0,0]
                Bwmm3_y = magnetic_field_p_Matrix_wmm[1,0]
                Bwmm3_z = magnetic_field_p_Matrix_wmm[2,0]
                
                Bwmm3_t = math.sqrt(Bwmm3_x**2 + Bwmm3_y**2  + Bwmm3_z**2)
               
               
                
                Biss_API_x = B_x-Bwmm3_x
                Biss_API_y = B_y-Bwmm3_y
                Biss_API_z = B_z-Bwmm3_z
                Biss_API_t = math.sqrt(Biss_API_x**2 + Biss_API_y**2  + Biss_API_z**2)
                
                
                # VEC Translation for verifing
                # Reverse Euler matrix 
                yawMatrix_API_ISS = np.matrix([
                [math.cos(1.57), -math.sin(1.57), 0],
                [math.sin(1.57), math.cos(1.57), 0],
                [0, 0, 1]
                ])

                pitchMatrix_API_ISS = np.matrix([
                [math.cos(pitch_api), 0, math.sin(pitch_api)],
                [0, 1, 0],
                [-math.sin(pitch_api), 0, math.cos(pitch_api)]
                ])

                rollMatrix_API_ISS = np.matrix([
                [1, 0, 0],
                [0, math.cos(roll_api), -math.sin(roll_api)],
                [0, math.sin(roll_api), math.cos(roll_api)]
                ])
                
                R_API_ISS = yawMatrix_API_ISS * pitchMatrix_API_ISS *  rollMatrix_API_ISS
                 
  
            
                # Magnetic field : processing 
                magnetic_field_p_Matrix = R_API_ISS *  magnetic_field_Matrix
                
                B_p_x = magnetic_field_p_Matrix[0,0]
                B_p_y = magnetic_field_p_Matrix[1,0]
                B_p_z = magnetic_field_p_Matrix[2,0]
                
                B_p_t = math.sqrt(B_p_x**2 + B_p_y**2  + B_p_z**2)
               
       
                
                # Direct Euler matrix WMM   NEC ---> VEC 
                pitch_iss = 0 
                roll_iss =  0 
                yaw_iss =  beta 
                
                yawMatrix_NEC_VEC= np.matrix([
                [math.cos(yaw_iss), math.sin(yaw_iss), 0],
                [-math.sin(yaw_iss), math.cos(yaw_iss), 0],
                [0, 0, 1]
                ])

                pitchMatrix_NEC_VEC= np.matrix([
                [math.cos(pitch_iss), 0, -math.sin(pitch_iss)],
                [0, 1, 0],
                [math.sin(pitch_iss), 0, math.cos(pitch_iss)]
                ])

                rollMatrix_NEC_VEC = np.matrix([
                [1, 0, 0],
                [0, math.cos(roll_iss), math.sin(roll_iss)],
                [0, -math.sin(roll_iss), math.cos(roll_iss)]
                ])
                
                R_NEC_VEC = rollMatrix_NEC_VEC * pitchMatrix_NEC_VEC *yawMatrix_NEC_VEC
                 
  
            
                # Magnetic field : processing
                magnetic_field_p_Matrix_wmm = R_NEC_VEC *  magnetic_field_Matrix_wmm
                
                Bwmm_p_x = magnetic_field_p_Matrix_wmm[0,0]
                Bwmm_p_y = magnetic_field_p_Matrix_wmm[1,0]
                Bwmm_p_z = magnetic_field_p_Matrix_wmm[2,0]
                
                Bwmm_p_t = math.sqrt(Bwmm_p_x**2 + Bwmm_p_y**2  + Bwmm_p_z**2)
               
               
                
                Biss_VEC_x = B_p_x-Bwmm_p_x
                Biss_VEC_y = B_p_y-Bwmm_p_y
                Biss_VEC_z = B_p_z-Bwmm_p_z
                Biss_VEC_t = math.sqrt(Biss_VEC_x**2 + Biss_VEC_y**2  + Biss_VEC_z**2)
                
                
                # FORCE - TORQUE - in NEC     
                
                """Coil Characteristics"""

                S= 1e-2
                N= 1.6 * 1e6
                Rez= 1             # superconductivity, for normal wire the resistence is about 3 MOhm
                l_bobina = 0.1
                I_sursa = 100     
                
                
                alpha = 0 # the angle between the coil and the horizontal plane
                
                gamma_x = (math.sin(alpha/2))**2 + math.cos(alpha)*(math.sin(beta/2))**2
                gamma_x = math.sqrt( gamma_x )
                gamma_x = 2 * math.asin(gamma_x)
                
                
                gamma_y = (math.sin(alpha/2))**2 + math.cos(alpha)*(math.sin((1.57-beta)/2))**2
                gamma_y = math.sqrt( gamma_y )
                gamma_y = 2 * math.asin(gamma_y)
                
                
                Bwmm_paralel = Bwmm_x * math.cos( gamma_x ) + Bwmm_y * math.cos( gamma_y )-Bwmm_z * math.sin( alpha) 
                Bwmm_perp = math.sqrt( Bwmm_t ** 2 - Bwmm_paralel ** 2)
                
                
                if line_count == 1: 
                    Bwmm_paralel_ant = Bwmm_paralel
                    
                vt_Bwmm_paralel =  (Bwmm_paralel - Bwmm_paralel_ant )/10
                
                Bwmm_paralel_ant = Bwmm_paralel
                
                dr = l_bobina * math.sin ( alpha )/1000
             
                lat_dr = lat + (l_bobina * math.cos ( alpha )* math.cos ( beta )) / r
                long_dr=long + (l_bobina * math.cos ( alpha )* math.sin ( beta )) / r
                
                 # Magetic field : WMM model  
                mag = gm.GeoMag(lat_dr,long_dr,408+dr, work_time)
                Bwmm_x2 = mag.bx / 1000
                Bwmm_y2 = mag.by / 1000
                Bwmm_z2 = mag.bz / 1000  
                Bwmm_t2 = math.sqrt(Bwmm_x2**2 + Bwmm_y2**2  + Bwmm_z2**2) 
                                
                 
                Bwmm_paralel2 = Bwmm_x2 * math.cos( gamma_x ) + Bwmm_y2 * math.cos( gamma_y )-Bwmm_z2 * math.sin( alpha) 
                
                vs_Bwmm_paralel = ( Bwmm_paralel2 - Bwmm_paralel)/l_bobina
                
                I_indus = - vt_Bwmm_paralel * S * N/ Rez *  10**(-6)
               
                
                mom_mag_indus = I_indus * S * N 
                mom_mag_sursa = I_sursa * S * N 
                
                Forta_wmm_indus = mom_mag_indus * vs_Bwmm_paralel * 10**(-6)
                Forta_wmm_sursa = mom_mag_sursa * vs_Bwmm_paralel * 10**(-6)
                
                
                mom_F_wmm_indus = mom_mag_indus *  Bwmm_perp * 10**(-6)
                mom_F_wmm_sursa = mom_mag_sursa *  Bwmm_perp * 10**(-6)
                
               
                # For the torque we calculated only the modulus, not its orientation.
                

               

                # print (rx, ry, rz)
                print("Row:  ", line_count)
              
                """
                print ("roll = ", roll)
                print ("pitch = ", pitch)
                print ("yaw = ", yaw)
                print ("B_achizitie :  ",B_x, B_y, B_z,B_t )
                print ("B_prelucrate : ",Bp_x, Bp_y, Bp_z, Bp_t)
                print ( "timp calculat:" , datetime_calc)
                print ( "timp curent=" ,  now_time)
                print ("lat = ", lat)
                print ("long = ", long)
                print ("DATA :  ",yr, mnh, dy, hr, mnt, sec, d_mjd2000 )
                
                print ("")
                print ("")
                """
                
               
                data_work= [ line_count, datetime_calc, pitch_rec, roll_rec, yaw_rec, beta, lat_rec_r, long_rec_r, lat_r, long_r, Theta_rad, B_x, B_y, B_z, B_t, Bth_x, Bth_y , Bth_z, Bth_t,Bwmm_x, Bwmm_y , Bwmm_z, Bwmm_t,Bwmm1_x, Bwmm1_y , Bwmm1_z,Bwmm2_x, Bwmm2_y , Bwmm2_z, Bwmm3_x, Bwmm3_y , Bwmm3_z, Bwmm3_t, pitch_api,roll_api, yaw_api, Biss_API_x, Biss_API_y , Biss_API_z, Biss_API_t, B_p_x, B_p_y, B_p_z, B_p_t, pitch_iss,roll_iss, yaw_iss,Bwmm_p_x, Bwmm_p_y , Bwmm_p_z, Bwmm_p_t, Biss_VEC_x, Biss_VEC_y , Biss_VEC_z, Biss_VEC_t, alpha, vt_Bwmm_paralel, vs_Bwmm_paralel, I_indus, mom_mag_indus, Forta_wmm_indus, mom_F_wmm_indus,  mom_mag_sursa, Forta_wmm_sursa, mom_F_wmm_sursa]
                data  = data_work
                csv_writer_f_output.writerow(data)
                line_count += 1






