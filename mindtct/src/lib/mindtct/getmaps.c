#include <stdio.h>
#include <lfs.h>
#include <mytime.h>
#include <log.h>

int lfs_detect_minutiae_V2_getmaps(int **odmap, int **olcmap, int **olfmap, int **ohcmap,
                        int *omw, int *omh,
                        unsigned char **obdata, int *obw, int *obh,
                        unsigned char *idata, const int iw, const int ih,
                        const LFSPARMS *lfsparms, int binarize)
{
   unsigned char *pdata, *bdata;
   int pw, ph, bw, bh;
   DIR2RAD *dir2rad;
   DFTWAVES *dftwaves;
   ROTGRIDS *dftgrids;
   ROTGRIDS *dirbingrids;
   int *direction_map, *low_contrast_map, *low_flow_map, *high_curve_map;
   int mw, mh;
   int ret, maxpad;

   /******************/
   /* INITIALIZATION */
   /******************/

   /* Determine the maximum amount of image padding required to support */
   /* LFS processes.                                                    */
   maxpad = get_max_padding_V2(lfsparms->windowsize, lfsparms->windowoffset,
                          lfsparms->dirbin_grid_w, lfsparms->dirbin_grid_h);

   /* Initialize lookup table for converting integer directions */
   /* to angles in radians.                                     */
   if((ret = init_dir2rad(&dir2rad, lfsparms->num_directions))){
      /* Free memory allocated to this point. */
      return(ret);
   }

   /* Initialize wave form lookup tables for DFT analyses. */
   /* used for direction binarization.                             */
   if((ret = init_dftwaves(&dftwaves, dft_coefs, lfsparms->num_dft_waves,
                        lfsparms->windowsize))){
      /* Free memory allocated to this point. */
      free_dir2rad(dir2rad);
      return(ret);
   }

   /* Initialize lookup table for pixel offsets to rotated grids */
   /* used for DFT analyses.                                     */
   if((ret = init_rotgrids(&dftgrids, iw, ih, maxpad,
                        lfsparms->start_dir_angle, lfsparms->num_directions,
                        lfsparms->windowsize, lfsparms->windowsize,
                        RELATIVE2ORIGIN))){
      /* Free memory allocated to this point. */
      free_dir2rad(dir2rad);
      free_dftwaves(dftwaves);
      return(ret);
   }

   /* Pad input image based on max padding. */
   if(maxpad > 0){   /* May not need to pad at all */
      if((ret = pad_uchar_image(&pdata, &pw, &ph, idata, iw, ih,
                             maxpad, lfsparms->pad_value))){
         /* Free memory allocated to this point. */
         free_dir2rad(dir2rad);
         free_dftwaves(dftwaves);
         free_rotgrids(dftgrids);
         return(ret);
      }
   }
   else{
      /* If padding is unnecessary, then copy the input image. */
      pdata = (unsigned char *)malloc(iw*ih);
      if(pdata == (unsigned char *)NULL){
         /* Free memory allocated to this point. */
         free_dir2rad(dir2rad);
         free_dftwaves(dftwaves);
         free_rotgrids(dftgrids);
         fprintf(stderr, "ERROR : lfs_detect_minutiae_V2 : malloc : pdata\n");
         return(-580);
      }
      memcpy(pdata, idata, iw*ih);
      pw = iw;
      ph = ih;
   }

   /* Scale input image to 6 bits [0..63] */
   /* !!! Would like to remove this dependency eventualy !!!     */
   /* But, the DFT computations will need to be changed, and     */
   /* could not get this work upon first attempt. Also, if not   */
   /* careful, I think accumulated power magnitudes may overflow */
   /* doubles.                                                   */
   bits_8to6(pdata, pw, ph);

   /******************/
   /*      MAPS      */
   /******************/

   /* Generate block maps from the input image. */
   if((ret = gen_image_maps(&direction_map, &low_contrast_map,
                    &low_flow_map, &high_curve_map, &mw, &mh,
                    pdata, pw, ph, dir2rad, dftwaves, dftgrids, lfsparms))){
      /* Free memory allocated to this point. */
      free_dir2rad(dir2rad);
      free_dftwaves(dftwaves);
      free_rotgrids(dftgrids);
      free(pdata);
      return(ret);
   }

   /* Deallocate working memories. */
   free_dir2rad(dir2rad);
   free_dftwaves(dftwaves);
   free_rotgrids(dftgrids);

   *odmap = direction_map;
   *olcmap = low_contrast_map;
   *olfmap = low_flow_map;
   *ohcmap = high_curve_map;
   *omw = mw;
   *omh = mh;

   if (!binarize) {
	   free(pdata);
	   return(0);
   }

   /******************/
   /* BINARIZARION   */
   /******************/

   /* Initialize lookup table for pixel offsets to rotated grids */
   /* used for directional binarization.                         */
   if((ret = init_rotgrids(&dirbingrids, iw, ih, maxpad,
                        lfsparms->start_dir_angle, lfsparms->num_directions,
                        lfsparms->dirbin_grid_w, lfsparms->dirbin_grid_h,
                        RELATIVE2CENTER))){
      /* Free memory allocated to this point. */
      free(pdata);
      free(direction_map);
      free(low_contrast_map);
      free(low_flow_map);
      free(high_curve_map);
      return(ret);
   }

   /* Binarize input image based on NMAP information. */
   if((ret = binarize_V2(&bdata, &bw, &bh,
                      pdata, pw, ph, direction_map, mw, mh,
                      dirbingrids, lfsparms))){
      /* Free memory allocated to this point. */
      free(pdata);
      free(direction_map);
      free(low_contrast_map);
      free(low_flow_map);
      free(high_curve_map);
      free_rotgrids(dirbingrids);
      return(ret);
   }

   /* Deallocate working memory. */
   free(pdata);
   free_rotgrids(dirbingrids);

   gray2bin(1, 1, 0, bdata, iw, ih);

   *obdata = bdata;
   *obw = bw;
   *obh = bh;

   return(0);
}

int get_maps(int **oquality_map,
                 int **odirection_map, int **olow_contrast_map,
                 int **olow_flow_map, int **ohigh_curve_map,
                 int *omap_w, int *omap_h,
                 unsigned char **obdata, int *obw, int *obh,
                 unsigned char *idata, const int iw, const int ih,
                 const int id, const double ppmm, const LFSPARMS *lfsparms, int binarize)
{
   int ret;
   int *direction_map, *low_contrast_map, *low_flow_map;
   int *high_curve_map, *quality_map;
   int map_w, map_h;
   unsigned char *bdata;
   int bw, bh;

   /* If input image is not 8-bit grayscale ... */
   if(id != 8){
      fprintf(stderr, "ERROR : get_minutiae : input image pixel ");
      fprintf(stderr, "depth = %d != 8.\n", id);
      return(-2);
   }

   /* Detect minutiae in grayscale fingerpeint image. */
   if((ret = lfs_detect_minutiae_V2_getmaps(&direction_map, &low_contrast_map,
                                   &low_flow_map, &high_curve_map,
                                   &map_w, &map_h,
                                   &bdata, &bw, &bh,
                                   idata, iw, ih, lfsparms, binarize))){
      return(ret);
   }

   /* Build integrated quality map. */
   if((ret = gen_quality_map(&quality_map,
                            direction_map, low_contrast_map,
                            low_flow_map, high_curve_map, map_w, map_h))){
      free(direction_map);
      free(low_contrast_map);
      free(low_flow_map);
      free(high_curve_map);
      return(ret);
   }

   /* Set output pointers. */
   *oquality_map = quality_map;
   *odirection_map = direction_map;
   *olow_contrast_map = low_contrast_map;
   *olow_flow_map = low_flow_map;
   *ohigh_curve_map = high_curve_map;
   *omap_w = map_w;
   *omap_h = map_h;
   *obdata = bdata;
   *obw = bw;
   *obh = bh;

   /* Return normally. */
   return(0);
}
