% eegplugin_std_dipoleDensity()
%
% Author: Makoto Miyakoshi JSPS/SCCN,INC,UCSD; Cincinnati Children's Hospital
%
% See also: eegplugin_std_dipoleDensity() std_pop_dipoleDensity() dipplot()        
%
% History
% 08/23/2024 Makoto and Komal. Re-visiting this plugin to make it work again.
% 03/06/2017 Makoto. Updated.
% 01/16/2015 ver 0.20 by Makoto. Changed function name.
% 02/19/2013 ver 1.1 by Makoto. STUDY = std_pop_dipoleDensity()
% 11/23/2012 ver 1.0 by Makoto. Created.

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2012, Makoto Miyakoshi JSPS/SCCN,INC,UCSD; Cincinnati Children's Hospital
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function eegplugin_std_dipoleDensity( fig, try_strings, catch_strings);

vers = 'std_dipoleDensity1.0';
if nargin < 3
    error('eegplugin_dipfit requires 3 arguments');
end

% create menu
std = findobj(fig, 'tag', 'study');
uimenu( std, 'label', 'Show cluster dipole densities', 'callback', 'STUDY = std_pop_dipoleDensity(STUDY, ALLEEG);', 'userdata', 'startup:off;study:on');