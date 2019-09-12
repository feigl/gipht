function M0inNm = okada_moment(shear_modulus,fault_length,fault_width,slip)
% calculate seismic moment M0 in Newton-meters
% 20190813 Kurt Feigl

M0inNm = shear_modulus * fault_length * fault_width * slip;
end

