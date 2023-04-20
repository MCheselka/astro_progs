#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#define BSZ			2
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

void main (int argc, char **argv)
{
	int		n, l, i, j, k, ptr, cnt;
	int		x, y, z, skip, sx, sy, fa;
	int		x1, y1, d1, d2, cx, cy;
	int		xxx, rrr, ggg, bbb;
	unsigned char	*uval1;
	float			*uval2, *uval3;
	float	min, max, sum1, sum2, utc;
	char	data[1024];
	FILE	*fd;
	color	*output;
	int		*mask, *idx;
	int		xmi, xma, ymi, yma;

    x = atoi (argv[3]);				// x size
    y = atoi (argv[4]);				// y size
    z = atoi (argv[5]);				// z size
    cx = atoi (argv[6]);			// start x location
    cy = atoi (argv[7]);			// start y location
    skip = atoi (argv[8]);			// image start number

	xmi = EDGE;
	xma = x - EDGE;
	ymi = EDGE;
	yma = y - EDGE;

	uval1 = (unsigned char *) calloc (x * y, sizeof (unsigned char));
	uval2 = (float *) calloc (x * y, sizeof (float));
	uval3 = (float *) calloc (x * y, sizeof (float));
	output = (color *) calloc (x * y, sizeof (color));

	if (cx < (x / 2)) {
		sx = cx + X_OFFSET;
	} else {
		sx = cx - X_OFFSET;
	}
	if (cy < (y / 2)) {
		sy = cy + Y_OFFSET;
	} else {
		sy = cy - Y_OFFSET;
	}

	mask = (int *) calloc (x * y, sizeof (int));

	srand (0);

	fd = fopen (argv[1], "rt");		// file list
	utc = 0;
	cnt = 0;
	l = 0;
	for (n = 0 ; n < z ; n++) {
		fscanf (fd, "%s", data);
		// fprintf (stderr,"%05d %s ",n,data);
		if (n >= skip) {
			read_fits (data, uval1, x, y, &utc);
			if ((cx != -1) && (cy != -1)) {
				// find peak
				max = (float) *(uval1 + ((cy * x) + cx));
				x1 = cx;
				y1 = cy;
				for (j = cy - BX_SZ ; j <= cy + BX_SZ ; j++) {
					for (i = cx - BX_SZ ; i <= cx + BX_SZ ; i++) {
						ptr = (j * x) + i;
						*(uval2+ptr) = (float) *(uval1+ptr);
						if (*(uval2+ptr) > max) {
							max = *(uval2+ptr);
							x1 = i;
							y1 = j;
						}
					}
				}
				// if it falls of the edge ...
				if ((x1 < xmi)||(x1 >= xma)||(y1 < ymi)||(y1 >= yma)||(sx < xmi)||(sx >= xma)||(sy < ymi)||(sy >= yma)) {
					break;
				}
				d1 = x1 - cx;
				d2 = y1 - cy;
				cx = x1;
				cy = y1;

				// do square aperture photometry on target
				sum1 = 0;
				for (j = cy - PH_SZ ; j <= cy + PH_SZ ; j++) {
					for (i = cx - PH_SZ ; i <= cx + PH_SZ ; i++) {
						ptr = (j * x) + i;
						*(uval2+ptr) = (float) *(uval1+ptr);
						*(uval3+ptr) += (float) *(uval1+ptr);
						sum1 += *(uval2+ptr);
						mask[ptr] += 1;
					}
				}

				// do square aperture photometry on sky (or another target)
				sum2 = 0;
				for (j = sy - PH_SZ ; j <= sy + PH_SZ ; j++) {
					for (i = sx - PH_SZ ; i <= sx + PH_SZ ; i++) {
						ptr = (j * x) + i;
						*(uval2+ptr) = (float) *(uval1+ptr);
						sum2 += *(uval2+ptr);
						// fprintf (stderr,"%d %d %f %f\n",i,j,*(uval2+ptr),sum2);
					}
				}
				sx += d1;
				sy += d2;

				// print out this stuff...
				// fprintf (stderr,"%8.4f %4d %4d %4d %4d %6.0f %6.0f %6.0f %03d %03d %03d\n",utc,cx[k],cy[k],d1,d2,sum1,sum2,sum1-sum2,rrr,ggg,bbb);

				fprintf (stderr,"%8.4f %4d %4d %4d %4d %6.0f %6.0f %6.0f\n",utc,cx,cy,d1,d2,sum1,sum2,sum1-sum2);

				// if the drift is too much, get a new object
				min = sqrt ((d1*d1)+(d2*d2));
				if (min > MIN_DIST) {
					break;
					// fprintf (stderr,"3 %d %d %f %d\n",n,k,utc,cnt);
					// fprintf (stderr,"#######################\n");
				}
			}
		}
	}
	fclose (fd);

	for (i = 0 ; i < x * y ; i++) {
		if (mask[i] > 0) {
			sum1 = uval3[i] / (float) mask[i];
			output[i].rr = (unsigned char) (sum1);
			output[i].gg = (unsigned char) (sum1);
			output[i].bb = (unsigned char) (sum1);
		}
	}

	fa = open (argv[2], O_RDWR | O_CREAT, S_IRWXU);
	write (fa, output, x * y * sizeof (color));
	close (fa);

	free (mask);
	free (output);
	free (uval1);
	free (uval2);
	free (uval3);
}
