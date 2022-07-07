### Tests
- [X] Test for all available data within 30 days took ~28 minutes (first 30 days of 2016).

### Todo
- [X] Check that modified codes give same results with Matteo's codes 
- [X] Check that our codes give similar results with Matlab codes 
- [X] Write code for manually picking Moho depths from migrated cross-sections. Modify so 
picks are made by pressing a button and not by clicking. Assing two types of
picks (i.e., one for Moho and one for uncertain moho)
- [X] Finish installation instrunctions
- [X] Read EPcrust velocity model. Make sure there are always Velocity values for sediments...
Checked that there is always vp and vp sedimentary.
- [X] Double-check back-azimuth changes along the theoretical wave path (fixed issue with using a function). 
- [ ] Compare results using EPcrust vs iasp91: 1) cross-sections 2) piercing points
- [X] Plot cross-section (update the colormap scale).
- [ ] Remove noisy RFs (see plot with individual trace rms values and choose a limit: 0.08?); Include this figure in the SI.
- [ ] Start ray paths from station elevation NOT z = 0 or z = -5...
- [ ] Double-check how swath is being calculated when creating 2D cross-sections!!!



### Secondary Todo
- [ ] Sort out codes for making RF plots.
- [ ] Organise codes `compute_RF_migration_spher.py`
- [ ] Finish tutorial
- [ ] Calculate TRF components...
