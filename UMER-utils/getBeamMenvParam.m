function [emittance, perveance, I, x0, y0, xp0, yp0, x0aper, y0aper, xp0aper, yp0aper] = getBeamMenvParam(BeamName)

switch BeamName
    case 'Pencil'
        %% Pencil Beam parameters
        emittance = 7.6; % mm-mrad (emit_x, emit_y)
        I = -0.6e-3; %ibeam current [Amps]
%         perveance = -4.9514e-09; %perveance calculated below
        % -- initial parameters
        x0        = 0.2; %cm
        y0        = 0.2; %cm
        xp0       = 0; %cm
        yp0       = 0; %cm
        % -- beam parameters at aperture
        x0aper      =  0.025;   %cm;
        y0aper      =  0.025;   %cm;
        xp0aper     = -0.00127; %cm
        yp0aper     = -0.00127; %cm
        
    case '6mA'
        %% 6mA Beam parameters
        emittance = 25.5; % mm-mrad (emit_x, emit_y)
        I = -6.0e-3; %ibeam current [Amps]
%         perveance = -4.9514e-08; %perveance calculated below
        % -- initial parameters
        x0        = 0.2; %cm
        y0        = 0.2; %cm
        xp0       = 0; %cm
        yp0       = 0; %cm
        % -- beam parameters at aperture
        x0aper      =  0.0875; %cm;
        y0aper      =  0.0875; %cm;
        xp0aper     = -0.0043; %??
        yp0aper     = -0.0043; %??
        
    case '23mA'
        %% 23mA Beam parameters
        emittance = 30; % mm-mrad (emit_x, emit_y)
        I = -21.0e-3; %ibeam current [Amps]
%         perveance =  -1.7330e-07; %perveance calculated below
        % -- initial parameters
        x0        = 0.2; %cm
        y0        = 0.2; %cm
        xp0       = 0; %cm
        yp0       = 0; %cm
        % -- beam parameters at aperture
        x0aper      =  0.15;    %cm;
        y0aper      =  0.15;    %cm;
        xp0aper     = -0.0067;  %cm
        yp0aper     = -0.0067;  %cm
        
    case '80mA'
        %% 80mA Beam parameters
        emittance = 58.9; % mm-mrad (emit_x, emit_y)
        I = -78.0e-3; %ibeam current [Amps]
%         perveance =  -6.4369e-07; %perveance calculated below
        % -- initial parameters
        x0        = 0.2; %cm
        y0        = 0.2; %cm
        xp0       = 0; %cm
        yp0       = 0; %cm
        % -- beam parameters at aperture
        x0aper      =  0.285;   %cm;
        y0aper      =  0.285;   %cm;
        xp0aper     = -0.01273; %cm
        yp0aper     = -0.01273; %cm

    case '100mA'
        %% 100mA Beam parameters
        emittance = 64.0; % mm-mrad (emit_x, emit_y)
        I = -104.0e-3; %ibeam current [Amps]
%         perveance = -8.5825e-07;
        %perveance calculated below
        % -- initial parameters
        x0        = 0.2; %cm
        y0        = 0.2; %cm
        xp0       = 0; %cm
        yp0       = 0; %cm
        % -- beam parameters at aperture
        x0aper      =  0.32;    %cm;
        y0aper      =  0.32;    %cm;
        xp0aper     = -0.0143;  %cm
        yp0aper     = -0.0143;  %cm

    otherwise
        error('User has not input a correct BeamName')
end
%% perveance Calculations
% From google
% m = 510.998946e3; %electron mass (eV/c^2)
% From Kiersten
% mo = 931.494e6; %atomic masss unit eV/c^2
Io = 17e3; %Budker (or Alfven) current in Amps. = 4pi(epsilon0)*((M*(c^3))/e)
% From Kappa2Current()
q = 1.602e-19;
M = 9.109e-31;
c = 2.9979e8;
E = 10.0e3;
T = q*E;
gamma = 1+T/(M*c^2);
beta = sqrt(1-gamma^-2);
% perveance = (q*lambda)/(2*pi()*epsilon0*M*(gamma^3)*(beta^2)*(c^2));%From slide sent by Kiersten
perveance = abs(2*((I)/Io)*(1/((gamma*beta)^3)));  %Simplification of above  
end

 
% """
%    Ubeams.py
%     Rami Kishek
%     Created: 8/2009
% 
%     Last Updated: 02/10/2011
% 
%   Beam parameters for UMER Beams
% 
%   Based on measurements by D. Stratakis, table 4.2 in PhD thesis
%   Current updated based on measurements Summer 2009 by S. Bernal
%   emittance for 23 mA beam updated based on measurements by H. Zhang
%   emittances for 80/100 mA beams extrapolated from that value
% 
%   operating points
%   I_Quad =
%        2.20, 2.068, 1.826, 1.518
% 
%   The matching is performed with nonlinear bgrd quads and dipoles
%   and using 40,000 particles
% 
%   TO DO:
%   !*! Matching needs to be redone using new magnet models (2010 version)
%   !*! preferably using TRACE to speed up the process
% """
% 
% UB = {
%     'pencil': {
%         'ibeam': -0.6e-3,
%         'emitx': 7.6e-6,
%         'emity': 7.6e-6,
%         'xcent': 0.0005,
%         'xpcent': 0.0,
%         'aper': [0.00025, 0.00025, -0.00127, -0.00127],
%         'ring': {'100%': [0.001650, 0.001371, -0.009274, 0.007927],
%                  '94%':  [0.001659, 0.001482, -0.008613, 0.007757],
%                  '83%':  [0.001717, 0.001663, -0.007678, 0.007346],
%                  '69%':  [0.001845, 0.001931, -0.006684, 0.006747],
%                 }
%         },
%     '6mA': {
%         'ibeam': -6.0e-3,
%         'emitx': 25.5e-6,
%         'emity': 25.5e-6,
%         'xcent': 0.0005,
%         'xpcent': 0.0,
%         'aper': [0.000875, 0.000875, -0.0043, -0.0043],
%         'ring': {'100%': [0.003218, 0.003085, -0.016961, 0.015836],
%                  '94%':  [0.003308, 0.003281, -0.016084, 0.015544],
%                  '83%':  [0.003533, 0.003679, -0.014735, 0.014916],
%                  '69%':  [0.003953, 0.004356, -0.013312, 0.014167],
%                 }
%         },
%     '23mA': {
%         'ibeam': -21.0e-3,
%         'emitx': 30.0e-6,
%         'emity': 30.0e-6,
%         'xcent': 0.0005,
%         'xpcent': 0.0,
%         'aper': [0.0015, 0.0015, -0.0067, -0.0067],
%         'ring': {'100%': [0.004578, 0.004775, -0.021763, 0.022386],
%                  '94%':  [0.004825, 0.005107, -0.021347, 0.022279],
%                  '83%':  [0.005364, 0.005803, -0.020576, 0.022015],
%                  '69%':  [0.006286, 0.007065, -0.019748, 0.021923],
%                 }
%         },
%     '80mA': {
%         'ibeam': -78.0e-3,
%         'emitx': 58.9e-6,
%         'emity': 58.9e-6,
%         'xcent': 0.0005,
%         'xpcent': 0.0,
%         'aper': [0.00285, 0.00285, -0.01273, -0.01273],
%         'ring': {'100%': [0.008366, 0.008814, -0.038605, 0.040407],
%                  '94%':  [0.008856, 0.009405, -0.038036, 0.040322],
%                  '83%':  [0.009916, 0.010764, -0.037323, 0.040331],
%                  '69%':  [0.011680, 0.013153, -0.036182, 0.040585],
%                 }
%         },
%     '100mA': {
%         'ibeam': -104.0e-3,
%         'emitx': 64.0e-6,
%         'emity': 64.0e-6,
%         'xcent': 0.0005,
%         'xpcent': 0.0,
%         'aper': [0.0032, 0.0032, -0.0143, -0.0143],
%         'ring': {'100%': [0.009570, 0.010084, -0.043915, 0.046068],
%                  '94%':  [0.010131, 0.010769, -0.043315, 0.046025],
%                  '83%':  [0.011369, 0.012309, -0.042482, 0.045980],
%                  '69%':  [0.013379, 0.015061, -0.041343, 0.046110],
%                 }
%         },
% }
% 
% class Beam(object):
%     """ class Beam(object): a container for beam parameters to pass them to input deck """
%     def __init__(self, extent, lEarth, aperture, oppoint):
%         b = UB[aperture]
%         self.ibeam = b['ibeam']
%         self.emitx = b['emitx']
%         self.emity = b['emity']
%         if lEarth and ('i' not in extent):
%             self.xcent  = b['xcent']
%             self.xpcent = b['xpcent']
%         else:
%             self.xcent  = 0.0
%             self.xpcent = 0.0
%         if 'i' in extent: ic = b['aper']
%         else:
%             ic = b['ring'][oppoint]
%         self.a0  = ic[0]
%         self.b0  = ic[1]
%         self.ap0 = ic[2]
%         self.bp0 = ic[3]
% 
%     def __repr__(self):
%         return """Loaded Beam parameters =
%     ibeam = %(ibeam)s
%     emitx = %(emitx)s
%     emity = %(emity)s
% 
%     xcent = %(xcent)s
%     xpcent = %(xpcent)s
% 
%     a0  = %(a0)s
%     ap0 = %(ap0)s
%     b0  = %(b0)s
%     bp0 = %(bp0)s
%     """ %vars(self)
