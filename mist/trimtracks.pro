;; this trims the MIST EEP tracks (not included) down to the essentials required by
;; EXOFASTv2 in order to load them faster and take less space. This
;; should only need to be run when updating the MIST tracks (and
;; should probably only be run by JDE -- if there are new MIST tracks
;; you'd like to use, please let me know.
;; NOTE: the tracks referenced here fix a few problems with the public
;; tracks as of 2018/02/07

pro trimtracks, v1=v1, mass_evolution=mass_evolution

;; prioritize files more likely to be used first
;;files = file_search('MIST_v1.1_tracks/feh_??.??*vvcrit0.?/eeps/?????M.track.eep', count=nfiles)
files = file_search(['MIST_v1.2_tracks/MIST_v1.2_feh_p?.??*vvcrit0.0_EEPS/?????M.track.eep',$
                     'MIST_v1.2_tracks/MIST_v1.2_feh_m?.??*vvcrit0.0_EEPS/?????M.track.eep',$
                     'MIST_v1.2_tracks/MIST_v1.2_feh_??.??*vvcrit0.4_EEPS/?????M.track.eep'],count=nfiles)

for i=0, nfiles-1 do begin

   filename=files[i] + '.idl'
   if file_test(filename) then continue
   
   eep_str = readmisteep(files[i])

   ;; log_surf_cell_z changed to log_surf_z at some point
   if (where(tag_names(eep_str) eq 'LOG_SURF_Z'))[0] eq -1 then logsurfz = eep_str.log_surf_cell_z $
   else logsurfz = eep_str.log_surf_z
   
   rstar = 10^eep_str.log_r                                       ;; r_sun
   teff = 10^eep_str.log_teff                                     ;; K
   age = eep_str.star_age/1d9                                     ;; Gyr
   feh = logsurfz - alog10(eep_str.surface_h1) - alog10(0.0181d0) ;; convert Z to [Fe/H], with Asplund 2009 (Z/H1)_sun=0.0181
   mass = eep_str.star_mass
   
   eep = dindgen(n_elements(rstar))+1d0
   ageweight = deriv(age,eep) ;; d(EEP)/d(Age) to transform the uniform EEP prior to a uniform Age prior 

   ;; structure better?
;   track = create_struct('age',age,'rstar',rstar,'teff',teff,'feh',feh,'mass',mass)
   track = transpose([[age],[rstar],[teff],[feh],[ageweight],[mass]])
   
   save, track, filename=filename

endfor

stop

end

