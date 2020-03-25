# ofdm_im
OFDM with Index Modulation (OFDM-IM) with various detection types and imperfect CSI
%% OFDM-IM simulation with ML, GD and LLR detectors and imperfect CSI
% Performance metics: SEP symbol error probability [1], BER in [2]
% Matlab version 2015b, also working well on 2019a.
% Note that MCIK-OFDM is another name of OFDM-IM see [1], [2].

%% Author information
% Thien Van Luong, Queen's University Belfast, UK, now with
% University of Southampton, UK.
% Email: tluong01@qub.ac.uk or thien.luong@soton.ac.uk.
% Personal page: https://tvluong.wordpress.com

%% References
% [1] T. V. Luong and Y. Ko, “A tight bound on BER of MCIK-OFDM with
% greedy detection and imperfect CSI,” IEEE Commun. Lett., vol. 21,
% no. 12, pp. 2594 – 2597, Dec. 2017.
% [2] T. V. Luong and Y. Ko, “Impact of CSI uncertainty on MCIK-OFDM:
% tight, closed-form symbol error probability analysis,” IEEE Trans. Veh.
% Technol., vol. 67, no. 2, pp. 1272 – 1279, Feb. 2018.

