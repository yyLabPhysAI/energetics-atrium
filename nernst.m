function voltage = nernst(concentration_in, concentration_out, valence)

voltage  = (RTONF/valence)*log(concentration_out/concentration_in);

end