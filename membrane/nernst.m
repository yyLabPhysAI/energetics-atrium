function voltage = nernst(concentration_in, concentration_out, valence, data)

voltage  = (data.RTONF/valence)*log(concentration_out/concentration_in);

end