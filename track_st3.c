#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#define BSZ			4
#define PSZ			9
#define BX_SZ		4
#define PH_SZ		3
#define EDGE		30
#define MAXNUM		10000
#define NUMOBJ		30
#define SIGMA		5
#define X_OFFSET	10
#define Y_OFFSET	10
#define MIN_DIST	5
#define MAXVAL		250

typedef struct {
	unsigned char	rr, gg, bb;
} color;

float norm_mean(float *val, int n)
{
    int x;
    float   s;

    s = 0;
    for (x = 0 ; x < n ; x++) {
        s += *(val+x);
        }

    return(s/n);
} /* norm_mean */

float standard_dev(float *val, int n)
{
    int x;
    double  s,t;

    s = t = 0;
    for (x = 0 ; x < n ; x++) {
        s += *(val+x) * *(val+x);
        t += *(val+x);
        }

    if (n == 1 || n == 0) {
        return((float)0);
    } else {
        return(sqrt(((s/(float)n)) - ((t/(float)n)*(t/(float)n))));
    }
}

void sort_it (float *data, int *idx, int num)
{
    int     i, j, m;

	// fprintf (stderr,"sorting %d items\n",num);
    // bubble sort over object 'size'
    for (j = 0 ; j < num ; j++) {
        for (i = 1 ; i < num - j ; i++) {
            if (data[idx[i-1]] <= data[idx[i]]) {
                m = idx[i];
                idx[i] = idx[i-1];
                idx[i-1] = m;
            }
        }
    }
}

float utc_dec (char *name)
{
	int		m;
    float   hh, mm, ss;
    float   utc_ret;
	char	*str, *token, nn[10];

	for (m = 0, str = name ; ; m++, str = NULL) {
		token = strtok (str, "_");
		if (token == NULL) {
			break;
		} else {
			// nothing
		}
		if (m == 3) {
			strncpy (nn+0,token+0,2);
			nn[2] = '\0';
			hh = atof (nn);
			strncpy (nn+0,token+2,2);
			nn[2] = '\0';
			mm = atof (nn);
			strncpy (nn+0,token+4,2);
			nn[2] = '\0';
			ss = atof (nn);
		}
	}

    utc_ret = (hh + (mm/60) + (ss/60/60));
    return (utc_ret);
}

int read_fits (char *name, unsigned char *data, int x, int y, float *utc)
{
    int     a, l, len, i;
    float   k;
    char    data1[82];
    char    data2[82];

    a = open (name, O_RDWR);

    memset (data, 0, 82);
    memset (data1, 0, 82);
    memset (data2, 0, 82);
    len = 0;
    while (read (a, data1, 80) != 0) {
        len++;
		// fprintf (stderr,"%s\n",data);
		if (strncmp (data1,"END",3) == 0) {
            break;
        }
		/*
        if (strncmp (data1,"UT",2) == 0) {
            strncpy (data2, data1+11, 8);
            // fprintf (stderr,"\tLST '%s'\n",data2);
            *utc = utc_dec (data2);
        }
		*/
        memset (data1, 0, 82);
        memset (data2, 0, 82);
    }
    *utc = utc_dec (name);

	// skip to the end of the header
	k = (float) len / 36.0;
    l = (int) (len / 36);
    if (k > (float) l) {
        l++;
    }
    l *= 2880;
	l += 2880;		// this is for GMN data
	// fprintf (stderr,"offset = %d\n",l);
	lseek (a, l, SEEK_SET);

    l = read (a, data, x * y * sizeof (unsigned char));

    close (a);
    return (l);
}

int get_st (unsigned char *data, int x, int y, int *cx, int *cy, int *sx, int *sy, int which, int numobj)
{
	int		i, j, k, n, m, ptr1, ptr2;
	int		tx, ty, cnt;
	int		*idx, *xx, *yy, *maski;
	float	*sums, *out;
	float	mean, stddev, thresh, max, sum;

	// fprintf (stderr,"get_st %05d %05d %04d %04d %04d %04d\n",which,numobj,cx[which],cy[which],sx[which],sy[which]);
	sums = (float *) calloc (MAXNUM, sizeof (float));
	xx = (int *) calloc (MAXNUM, sizeof (int));
	yy = (int *) calloc (MAXNUM, sizeof (int));
	idx = (int *) calloc (MAXNUM, sizeof (int));
	out = (float *) calloc (x * y, sizeof (float));
	maski = (int *) calloc (x * y, sizeof (int));

	for (i = 0 ; i < x * y ; i++) {
		*(out+i) = (float) *(data+i);
	}

	for (k = 0 ; k < MAXNUM ; k++) {
		idx[k] = k;
		xx[k] = -1;
		yy[k] = -1;
	}

	// when which == numobj, this is a new run and there are no masks
	if (which != numobj) {
		for (k = 0 ; k < numobj ; k++) {
			// mask off all objects in the current image that are going to be kept
			if ((k != which) && (cx[k] != -1) && (cy[k] != -1)) {
				// fprintf (stderr,"blocking %04d %04d %04d\n",k,cx[k],cy[k]);
				for (n = cy[k] - PSZ ; n <= cy[k] + PSZ ; n++) {
					for (m = cx[k] - PSZ ; m <= cx[k] + PSZ ; m++) {
						ptr2 = (n * x) + m;
						maski[ptr2] = 1;
					}
				}
			}
		}
	}

	// get some image stats (maybe do this a different way?)
	mean = norm_mean (out, x * y);
	stddev = standard_dev (out, x * y);
	thresh = mean + (stddev * SIGMA);

	// find peaks based on thresholding
	cnt = 0;
	// fprintf (stderr,"%s %f %f %f\n",data,utc,mean,stddev);
	for (j = EDGE ; j < (y - EDGE) ; j++) {
		for (i = EDGE ; i < (x - EDGE) ; i++) {
			ptr1 = (j * x) + i;
			if ((maski[ptr1] == 0) && (out[ptr1] > thresh) && (out[ptr1] < MAXVAL)) {
				max = out[ptr1];
				tx = i;
				ty = j;
				for (n = j - BSZ ; n <= j + BSZ ; n++) {
					for (m = i - BSZ ; m <= i + BSZ ; m++) {
						ptr2 = (n * x) + m;
						if (out[ptr2] > max) {
							max = out[ptr2];
							tx = m;
							ty = n;
						}
					}
				}
				sum = 0;
				for (n = ty - PSZ ; n <= ty + PSZ ; n++) {
					for (m = tx - PSZ ; m <= tx + PSZ ; m++) {
						ptr2 = (n * x) + m;
						sum += out[ptr2];
						maski[ptr2] = 1;
					}
				}
				// fprintf (stderr,"loc %04d %04d %04d %8.2f\n",cnt,tx,ty,sum);
				sums[cnt] = sum;
				xx[cnt] = tx;
				yy[cnt] = ty;
				cnt++;
			}
		}
	}

	sort_it (sums, idx, cnt);


	/*
	fprintf (stderr,"%f %d\n",utc,cnt);
	for (k = 0 ; k < cnt ; k++) {
		fprintf (stderr,"%04d %04d %.2f %04d %04d\n",k,idx[k],sums[idx[k]],xx[idx[k]],yy[idx[k]]);
	}
	*/

	if (which == numobj) {
		// assign them all
		for (k = 0 ; k < cnt ; k++) {
			cx[k] = xx[idx[k]];
			cy[k] = yy[idx[k]];
			if (xx[idx[k]] < x/2) {
				sx[k] = xx[idx[k]] + X_OFFSET;
			} else {
				sx[k] = xx[idx[k]] - X_OFFSET;
			}
			if (yy[idx[k]] < y/2) {
				sy[k] = yy[idx[k]] + Y_OFFSET;
			} else {
				sy[k] = yy[idx[k]] - Y_OFFSET;
			}
			// fprintf (stderr,"assign1 %04d %04d %.2f %04d %04d\n",k,idx[k],sums[idx[k]],xx[idx[k]],yy[idx[k]]);
		}
		for (k = cnt ; k < MAXNUM ; k++) {
			cx[k] = sx[k] = xx[idx[k]];
			cy[k] = sy[k] = yy[idx[k]];
		}
	} else {
		// assign just the brighest one
		cx[which] = xx[idx[0]];
		cy[which] = yy[idx[0]];
		if (cx[which] < (x/2)) {
			sx[which] = cx[which] + X_OFFSET;
		} else {
			sx[which] = cx[which] - X_OFFSET;
		}
		if (yy[idx[0]] < (y/2)) {
			sy[which] = cy[which] + Y_OFFSET;
		} else {
			sy[which] = cy[which] - Y_OFFSET;
		}
		// fprintf (stderr,"assign2 %04d %04d %.2f %04d %04d %04d %04d\n",which,idx[0],sums[idx[0]],cx[which],cy[which],sx[which],sy[which]);
	}

	free (maski);
	free (out);
	free (sums);
	free (xx);
	free (yy);
	free (idx);

	return (cnt);
}

void main (int argc, char **argv)
{
	int		n, l, i, j, k, ptr, cnt;
	int		x, y, z, numobj, *sx, *sy, fa;
	int		x1, y1, d1, d2, *cx, *cy;
	int		xxx, rrr, ggg, bbb;
	int		o1;
	unsigned char	*uval1;
	float			*uval2, *uval3;
	float	min, max, sum1, sum2, utc;
	char	data[1024];
	FILE	*fd, *ft;
	color	*output;
	int		*mask, *mask1, *idx;
	int		xmi, xma, ymi, yma;

    x = atoi (argv[2]);				// x size
    y = atoi (argv[3]);				// y size
    z = atoi (argv[4]);				// z size
    numobj = atoi (argv[5]);		// number of stars to track at one time
    o1 = atoi (argv[6]);			// which star to mark with green
	/*
	cx = atoi (argv[5]);			// search x center
	cy = atoi (argv[6]);			// search y center
	sx = atoi(argv[7]);				// sky x location
	sy = atoi(argv[8]);				// sky y location
	*/

	xmi = EDGE;
	xma = x - EDGE;
	ymi = EDGE;
	yma = y - EDGE;
	// fprintf (stderr,"%d %d %d %d\n",xmi,xma,ymi,yma);

	uval1 = (unsigned char *) calloc (x * y, sizeof (unsigned char));
	uval2 = (float *) calloc (x * y, sizeof (float));
	uval3 = (float *) calloc (x * y, sizeof (float));
	output = (color *) calloc (x * y, sizeof (color));

	cx = (int *) calloc (MAXNUM, sizeof (int));
	cy = (int *) calloc (MAXNUM, sizeof (int));
	sx = (int *) calloc (MAXNUM, sizeof (int));
	sy = (int *) calloc (MAXNUM, sizeof (int));
	mask = (int *) calloc (x * y, sizeof (int));
	mask1 = (int *) calloc (x * y, sizeof (int));

	srand (0);

	ft = fopen ("starcount", "wt");

	fd = fopen (argv[1], "rt");		// file list
	utc = 0;
	cnt = 0;
	l = 0;
	for (n = 0 ; n < z ; n++) {
		fscanf (fd, "%s", data);
		// fprintf (stderr,"%05d %s ",n,data);
		read_fits (data, uval1, x, y, &utc);
		if (n == 0) {
			// get the first set of objects to track
			cnt += get_st (uval1, x, y, cx, cy, sx, sy, numobj, numobj);
			l++;
			// fprintf (stderr,"1 %d %f %d\n",n,utc,cnt);
			/*
			for (k = 0 ; k < numobj ; k++) {
				fprintf (stderr,"--- %8.4f %4d %4d %4d %4d\n",utc,cx[k],cy[k],sx[k],sy[k]);
			}
			*/
		} else {
			for (k = 0 ; k < numobj ; k++) {
				if ((cx[k] != -1) && (cy[k] != -1)) {
					x1 = cx[k];
					y1 = cy[k];
					// if limits are reached on any object, find a new object
					if ((x1 < xmi) || (x1 >= xma) || (y1 < ymi) || (y1 >= yma) || (sx[k] < xmi) || (sx[k] >= xma) || (sy[k] < ymi) || (sy[k] >= yma)) {
						cnt += get_st (uval1, x, y, cx, cy, sx, sy, k, numobj);
						l++;
						// fprintf (stderr,"2 %d %d %f %d\n",n,k,utc,cnt);
						// fprintf (stderr,"#######################\n");
					}

					// find peak
					max = (float) *(uval1 + ((cy[k] * x) + cx[k]));
					x1 = cx[k];
					y1 = cy[k];
					for (j = cy[k] - BX_SZ ; j <= cy[k] + BX_SZ ; j++) {
						for (i = cx[k] - BX_SZ ; i <= cx[k] + BX_SZ ; i++) {
							ptr = (j * x) + i;
							*(uval2+ptr) = (float) *(uval1+ptr);
							if (*(uval2+ptr) > max) {
								max = *(uval2+ptr);
								x1 = i;
								y1 = j;
							}
						}
					}
					d1 = x1 - cx[k];
					d2 = y1 - cy[k];
					cx[k] = x1;
					cy[k] = y1;
					if ((k == o1) && (o1 >= 0)) {
						ptr = (cy[k] * x) + cx[k];
						mask1[ptr] = 1;
					}

					// do square aperture photometry on target
					sum1 = 0;
					for (j = cy[k] - PH_SZ ; j <= cy[k] + PH_SZ ; j++) {
						for (i = cx[k]- PH_SZ ; i <= cx[k] + PH_SZ ; i++) {
							ptr = (j * x) + i;
							*(uval2+ptr) = (float) *(uval1+ptr);
							*(uval3+ptr) += (float) *(uval1+ptr);
							sum1 += *(uval2+ptr);
							mask[ptr] += 1;
						}
					}

					// do square aperture photometry on sky (or another target)
					sum2 = 0;
					for (j = sy[k] - PH_SZ ; j <= sy[k] + PH_SZ ; j++) {
						for (i = sx[k] - PH_SZ ; i <= sx[k] + PH_SZ ; i++) {
							ptr = (j * x) + i;
							*(uval2+ptr) = (float) *(uval1+ptr);
							sum2 += *(uval2+ptr);
							// fprintf (stderr,"%d %d %f %f\n",i,j,*(uval2+ptr),sum2);
						}
					}
					sx[k] += d1;
					sy[k] += d2;

					// print out this stuff...
					// fprintf (stderr,"%8.4f %4d %4d %4d %4d %6.0f %6.0f %6.0f %03d %03d %03d\n",utc,cx[k],cy[k],d1,d2,sum1,sum2,sum1-sum2,rrr,ggg,bbb);

					// HERE
					// if ((k == o1) && (o1 >= 0)) {
						fprintf (stderr,"%8.4f %4d %4d %4d %4d %6.0f %6.0f %6.0f\n",utc,cx[k],cy[k],d1,d2,sum1,sum2,sum1-sum2);
					// }

					// if the drift is too much, get a new object
					min = sqrt ((d1*d1)+(d2*d2));
					if (min > MIN_DIST) {
						cnt += get_st (uval1, x, y, cx, cy, sx, sy, k, numobj);
						l++;
						// fprintf (stderr,"3 %d %d %f %d\n",n,k,utc,cnt);
						// fprintf (stderr,"#######################\n");
					}
	
					// fprintf (stderr,"%s %d %d %f %f\n",data,cx,cy,sum1,sum2);
					// fprintf (stderr,"%s %d %d %f\n",data,sx,sy,sum1-sum2);
				} else {
					// fprintf (stderr,"# %04d %04d\n",n,k);
					cnt += get_st (uval1, x, y, cx, cy, sx, sy, k, numobj);
					l++;
					// fprintf (stderr,"4 %d %d %f %d\n",n,k,utc,cnt);
				}
			}
		}
		if (l != 0) {
			fprintf (ft,"%f %f\n",utc,(float)cnt/(float)l);
		}
	}
	fclose (fd);
	fclose (ft);

	for (i = 0 ; i < x * y ; i++) {
		if (mask[i] > 0) {
			sum1 = uval3[i] / (float) mask[i];
			if (mask1[i] > 0) {
				output[i].rr = (unsigned char) 0;
				output[i].gg = (unsigned char) (sum1);
				output[i].bb = (unsigned char) 0;
			} else {
				output[i].rr = (unsigned char) (sum1);
				output[i].gg = (unsigned char) (sum1);
				output[i].bb = (unsigned char) (sum1);
			}
			/*
			if (sum1 <= 86) {
				output[i].rr = (unsigned char) sum1;
				output[i].gg = (unsigned char) 0;
				output[i].bb = (unsigned char) 0;
			} else if ((sum1 > 86) && (sum1 <= 172)) {
				output[i].rr = (unsigned char) 0;
				output[i].gg = (unsigned char) sum1;
				output[i].bb = (unsigned char) 0;
			} else if (sum1 > 172) {
				output[i].rr = (unsigned char) 0;
				output[i].gg = (unsigned char) 0;
				output[i].bb = (unsigned char) sum1;
			}
			*/

			/*
			j = (int)(((float)rand() / (float)RAND_MAX) * (uval3[i] / (float) mask[i]));
			output[i].rr = (unsigned char) j;
			j = (int)(((float)rand() / (float)RAND_MAX) * (uval3[i] / (float) mask[i]));
			output[i].gg = (unsigned char) j;
			j = (int)(((float)rand() / (float)RAND_MAX) * (uval3[i] / (float) mask[i]));
			output[i].bb = (unsigned char) j;
			*/
		}
	}

	fa = open ("long_exp2.data", O_RDWR | O_CREAT, S_IRWXU);
	write (fa, output, x * y * sizeof (color));
	close (fa);

	free (cx);
	free (cy);
	free (sx);
	free (sy);
	free (mask1);
	free (mask);
	free (output);
	free (uval1);
	free (uval2);
	free (uval3);
}
