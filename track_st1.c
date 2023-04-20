#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#define BX_SZ	4
#define PH_SZ	3
#define EDGE	10

typedef struct {
	unsigned char	rr, gg, bb;
} color;

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
	int		n, l, i, j, ptr;
	int		x, y, z, sx, sy, fa;
	int		x1, y1, d1, d2, c_x1, c_y1;
	unsigned char	*uval1;
	float			*uval2, *uval3;
	float	min, max, sum1, sum2, utc;
	char	data[1024];
	FILE	*fd;
	color	*output;
	int		*mask;
	int		xmi, xma, ymi, yma;

	fd = fopen (argv[1], "rt");		// file list
    x = atoi (argv[2]);				// x size
    y = atoi (argv[3]);				// y size
    z = atoi (argv[4]);				// z size
	c_x1 = atoi (argv[5]);			// search x center
	c_y1 = atoi (argv[6]);			// search y center
	sx = atoi(argv[7]);				// sky x location
	sy = atoi(argv[8]);				// sky y location

	xmi = EDGE;
	xma = x - EDGE;
	ymi = EDGE;
	yma = y - EDGE;
	fprintf (stderr,"%d %d %d %d\n",xmi,xma,ymi,yma);

	uval1 = (unsigned char *) calloc (x * y, sizeof (unsigned char));
	uval2 = (float *) calloc (x * y, sizeof (float));
	uval3 = (float *) calloc (x * y, sizeof (float));
	output = (color *) calloc (x * y, sizeof (color));

	mask = (int *) calloc (x * y, sizeof (int));
	for (i = 0 ; i < x * y ; i++) {
		mask[i] = 1;
	}

	utc = 0;
	for (n = 0 ; n < z ; n++) {
		x1 = c_x1;
		y1 = c_y1;
		if ((x1 < xmi) || (x1 >= xma) || (y1 < ymi) || (y1 >= yma)) {
			// fprintf (stderr,"1 %d %d %d %d %d %d\n",x1,y1,xmi,xma,ymi,yma);
			continue;
		}
		if ((sx < xmi) || (sx >= xma) || (sy < ymi) || (sy >= yma)) {
			// fprintf (stderr,"2 %d %d %d %d %d %d\n",sx,sy,xmi,xma,ymi,yma);
			continue;
		}
		fscanf (fd, "%s", data);
		l = read_fits (data, uval1, x, y, &utc);
		// find peak
		max = (float) *(uval1 + ((c_y1 * x) + c_x1));
		for (j = c_y1 - BX_SZ ; j <= c_y1 + BX_SZ ; j++) {
			for (i = c_x1 - BX_SZ ; i <= c_x1 + BX_SZ ; i++) {
				ptr = (j * x) + i;
				*(uval2+ptr) = (float) *(uval1+ptr);
				// fprintf (stderr,"%03d ",((int)uval1[ptr] % 8));
				// fprintf (stderr,"%03d %03d %03d\n",i,j,(int)uval1[ptr]);
				if (*(uval2+ptr) > max) {
					max = *(uval2+ptr);
					x1 = i;
					y1 = j;
				}
			}
			// fprintf (stderr,"\n");
		}
		d1 = x1 - c_x1;
		d2 = y1 - c_y1;
		/*
		d1 = c_x1 - x1;
		d2 = c_y1 - y1;
		*/
		c_x1 = x1;
		c_y1 = y1;

		// do square aperture photometry on target
		sum1 = 0;
		for (j = c_y1 - PH_SZ ; j <= c_y1 + PH_SZ ; j++) {
			for (i = c_x1- PH_SZ ; i <= c_x1 + PH_SZ ; i++) {
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
		fprintf (stderr,"%8.4f %4d %4d %4d %4d %6.0f %6.0f %6.0f\n",utc,c_x1,c_y1,d1,d2,sum1,sum2,sum1-sum2);
		// fprintf (stderr,"%s %d %d %f %f\n",data,c_x1,c_y1,sum1,sum2);
		// fprintf (stderr,"%s %d %d %f\n",data,sx,sy,sum1-sum2);
	}
	fclose (fd);

	for (i = 0 ; i < x * y ; i++) {
		if (mask[i] > 0) {
			output[i].rr = (unsigned char) (uval3[i] / (float) mask[i]);
			output[i].gg = (unsigned char) (uval3[i] / (float) mask[i]);
			output[i].bb = (unsigned char) (uval3[i] / (float) mask[i]);
		}
	}

	fa = open ("long_exp.data", O_RDWR | O_CREAT, S_IRWXU);
	write (fa, output, x * y * sizeof (color));
	close (fa);

	free (mask);
	free (output);
	free (uval1);
	free (uval2);
	free (uval3);
}
