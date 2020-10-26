function tdb2tcb, tdb

lb = 1.550519768d-8
tdb0 = -6.55d-5
t0 = 2443144.5003725d0

tcb = tdb + ((tdb - t0)*lb*86400d0 + tdb0)/86400.d0

return, tcb

end
