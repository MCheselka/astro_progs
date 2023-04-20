#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#define BX_SZ	2
#define PH_SZ	8
#define EDGE	10
#define SIGMA	5

typedef struct {
	unsigned char	rr, gg, bb;
} color;

void sort_it (float *data, int *idx, int num)
{
	int		i, j, m;

	/*
	for (i = 0 ; i < obj_num ; i++) {
		idx[i] = i;
	}
	*/
    
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

// int compare (float a, float b)
int compare (const void * a, const void * b)
{
	int		rv;

	rv = 0;

	if (*(float *)a > *(float *)b) {
		rv = 1;
	} else if (*(float *)a < *(float *)b) {
		rv = -1;
	}
	return (rv);
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

	/*
    sscanf (utc, "%f:%f:%f",&hh,&mm,&ss);

	*/

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
	int		m, n, l, i, j, k, ptr1, ptr2;
	int		x, y, z, sx, sy, fa, cnt;
	unsigned char	*uval1;
	float			*uval2, *uval3;
	float	mean, stddev, thresh, max, utc, sum, *sums;
	char	data[1024];
	FILE	*fd;
	color	*output;
	int		*mask, *idx, *xx, *yy;
	int		xmi, xma, ymi, yma;

	fd = fopen (argv[1], "rt");		// file list
    x = atoi (argv[2]);				// x size
    y = atoi (argv[3]);				// y size
    z = atoi (argv[4]);				// z size

	xmi = EDGE;
	xma = x - EDGE;
	ymi = EDGE;
	yma = y - EDGE;
	// fprintf (stderr,"%d %d %d %d\n",xmi,xma,ymi,yma);

	uval1 = (unsigned char *) calloc (x * y, sizeof (unsigned char));
	uval2 = (float *) calloc (x * y, sizeof (float));
	uval3 = (float *) calloc (x * y, sizeof (float));
	sums = (float *) calloc (x * y, sizeof (float));
	output = (color *) calloc (x * y, sizeof (color));

	mask = (int *) calloc (x * y, sizeof (int));
	idx = (int *) calloc (x * y, sizeof (int));
	xx = (int *) calloc (x * y, sizeof (int));
	yy = (int *) calloc (x * y, sizeof (int));

	utc = 0;
	update = 10;
	for (k = 0 ; k < z ; k++) {
		fscanf (fd, "%s", data);
		l = read_fits (data, uval1, x, y, &utc);
		for (i = 0 ; i < x * y ; i++) {
			*(uval2+i) = (float) *(uval1+i);
		}
		mean = norm_mean (uval2, x * y);
		stddev = standard_dev (uval2, x * y);
		thresh = mean + (stddev * SIGMA);
		cnt = 0;
		// fprintf (stderr,"%s %f %f %f\n",data,utc,mean,stddev);
		for (j = EDGE ; j < (y - EDGE) ; j++) {
			for (i = EDGE ; i < (x - EDGE) ; i++) {
				ptr1 = (j * x) + i;
				if (mask[ptr1] == 0) {
					if (uval2[ptr1] > thresh) {
						max = uval2[ptr1];
						sx = i;
						sy = j;
						for (n = j - BX_SZ ; n <= j + BX_SZ ; n++) {
							for (m = i - BX_SZ ; m <= i + BX_SZ ; m++) {
								ptr2 = (n * x) + m;
								if (uval2[ptr2] > max) {
									max = uval2[ptr2];
									sx = m;
									sy = n;
								}
							}
						}
						sum = 0;
						for (n = sy - PH_SZ ; n <= sy + PH_SZ ; n++) {
							for (m = sx - PH_SZ ; m <= sx + PH_SZ ; m++) {
								ptr2 = (n * x) + m;
								uval3[ptr2] = uval2[ptr2];
								sum += uval3[ptr2];
								mask[ptr2] = 1;
							}
						}
						// fprintf (stderr,"%04d %04d %04d %8.2f\n",cnt,sx,sy,sum);
						sums[cnt] = sum;
						xx[cnt] = sx;
						yy[cnt] = sy;
						cnt++;
					}
				}
			}
		}

		for (k = 0 ; k < cnt ; k++) {
			idx[k] = k;
			fprintf (stderr,"%04d %04d %.2f %04d %04d\n",k,idx[k],sums[idx[k]],xx[idx[k]],yy[idx[k]]);
		}
		sort_it (sums, idx, cnt);
		for (k = 0 ; k < cnt ; k++) {
			fprintf (stderr,"%04d %04d %.2f %04d %04d\n",k,idx[k],sums[idx[k]],xx[idx[k]],yy[idx[k]]);
		}
	}
	fclose (fd);

	for (i = 0 ; i < x * y ; i++) {
		output[i].rr = (unsigned char) (uval3[i] * mask[i]);
		output[i].gg = (unsigned char) (uval3[i] * mask[i]);
		output[i].bb = (unsigned char) (uval3[i] * mask[i]);
	}
	fa = open ("objects.data", O_RDWR | O_CREAT, S_IRWXU);
	write (fa, output, x * y * sizeof (color));
	close (fa);

	/*
	for (i = 0 ; i < x * y ; i++) {
		output[i].rr = (unsigned char) (uval3[i]);
		output[i].gg = (unsigned char) (uval3[i]);
		output[i].bb = (unsigned char) (uval3[i]);
	}
	fa = open ("objects.data", O_RDWR | O_CREAT, S_IRWXU);
	write (fa, output, x * y * sizeof (color));
	close (fa);
	*/

	free (sums);
	free (xx);
	free (yy);
	free (idx);
	free (mask);
	free (output);
	free (uval1);
	free (uval2);
	free (uval3);
}
