#include "ReShadeUI.fxh"

uniform int iKuwahara_Radius < __UNIFORM_SLIDER_INT1
	ui_min = 1; ui_max = 16;
	ui_label = "Radius";
	ui_tooltip = "Silence!";
> = 5;

#include "ReShade.fxh"

float fmod(float x, float y) {
    return x - y * floor(x/y);
}

void PS_Kuwahara_Basic(float4 vpos : SV_Position, float2 texcoord : TEXCOORD, out float4 res : SV_Target0)
{
    int radius = iKuwahara_Radius;
    float n = float((radius+1)*(radius+1));

    float3 means[4] = {float3(0,0,0), float3(0,0,0), float3(0,0,0), float3(0,0,0)};
    float3 sigma[4] = {float3(0,0,0), float3(0,0,0), float3(0,0,0), float3(0,0,0)};

    int4 bounds[4] = {
        int4(-radius, 0, 0, radius),
        int4(0, radius, 0, radius),
        int4(-radius, 0, -radius, 0),
        int4(0, radius, -radius, 0)
    };

    for (int k=0; k<4; k++) {
        for (int i=bounds[k].x; i<=bounds[k].y; i++) {
            for (int j=bounds[k].z; j<=bounds[k].w; j++) {
                float2 offset = BUFFER_PIXEL_SIZE*float2(i,j);

                float3 val = tex2D(ReShade::BackBuffer, texcoord + offset).rgb;
                
                //For now this contains sums
                means[k]  += val;
                //And sums of squares
                sigma[k] += val*val;
            }
        }
    }

    //Large number
    float min_sigma = 1000.0;
    float3 color = float3(0,0,0);
    
    for (int k=0; k<4; k++) {
        //Now it's actually means
        means[k] /= n;
        //And standard deviations
        sigma[k] = abs(sigma[k]/n - means[k]*means[k]);
        
        //Sum of standard deviations for each color channel
        float sigma2 = sigma[k].x + sigma[k].y + sigma[k].z;
        
        //Use mean from region with lowest deviation
        if (sigma2 < min_sigma) {
            min_sigma = sigma2;
            color = means[k];
        }
    }

    //=====Output:
	res.xyz = color;
	res.w = 1.0;
}

void PS_Kuwahara_Generalized(float4 vpos : SV_Position, float2 texcoord : TEXCOORD, out float4 res : SV_Target0) {
    const float PI  = 3.1415926535;
    const float GOLDEN_RATIO = 1.6180339887;

    float3 color = float3(0,0,0);

    //Number equivalent to standard technique
    const int samples = 4*(iKuwahara_Radius+1)*(iKuwahara_Radius+1);
    const float radius = float(iKuwahara_Radius);

    float3 means[8] = {float3(0,0,0), float3(0,0,0), float3(0,0,0), float3(0,0,0),
                       float3(0,0,0), float3(0,0,0), float3(0,0,0), float3(0,0,0)};

    float3 sigma[8] = {float3(0,0,0), float3(0,0,0), float3(0,0,0), float3(0,0,0),
                       float3(0,0,0), float3(0,0,0), float3(0,0,0), float3(0,0,0)};

    float norms[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    //Sample from disk, along a Fibbonaci spiral:
    const float2 scale = radius * BUFFER_PIXEL_SIZE;
    const float dt = 1.0/float(samples);
    const float A = 2.0*PI*float(samples)/GOLDEN_RATIO;

    float t = 0.0;

    for (int i=0; i<samples; i++) {
        //Generate point on the spiral:
        float ang = A*t;
        float2 offset = scale * sqrt(t) * float2(sin(ang), cos(ang));

        //Retrieve index of the quadrant
        float f_ID = fmod(4.0*ang/PI, 8.0);
        int i_ID = int(f_ID);

        //Gather values
        float3 val = tex2D(ReShade::BackBuffer, texcoord + offset).rgb;
                
        //For now this contains sums
        means[i_ID] += val;
        //And sums of squares
        sigma[i_ID] += val*val;
        //We now also need normalizations
        norms[i_ID] += 1.0;

        //Increment along the spiral
        t += dt;
    }

    //Large number
    float min_sigma = 1000.0;
    
    for (int k=0; k<8; k++) {
        //Now it's actually means
        means[k] /= norms[k];
        //And standard deviations
        sigma[k] = abs(sigma[k]/norms[k] - means[k]*means[k]);
        
        //Sum of standard deviations for each color channel
        float sigma2 = sigma[k].x + sigma[k].y + sigma[k].z;
        
        //Use mean from region with lowest deviation
        if (sigma2 < min_sigma) {
            min_sigma = sigma2;
            color = means[k];
        }
    }

    //=====Output:
	res.xyz = color;
	res.w = 1.0;
}

void PS_Kuwahara_Anisotropic(float4 vpos : SV_Position, float2 texcoord : TEXCOORD, out float4 res : SV_Target0) {
    float2 d = BUFFER_PIXEL_SIZE;

    //Approximation of partial derivatives by using the Sobel operator:
    float3 f_x = -1.0 * tex2D(ReShade::BackBuffer, texcoord + float2(-d.x, -d.y)).xyz
                 -2.0 * tex2D(ReShade::BackBuffer, texcoord + float2(-d.x,  0.0)).xyz
                 -1.0 * tex2D(ReShade::BackBuffer, texcoord + float2(-d.x,  d.y)).xyz
                 +1.0 * tex2D(ReShade::BackBuffer, texcoord + float2( d.x, -d.y)).xyz
                 +2.0 * tex2D(ReShade::BackBuffer, texcoord + float2( d.x,  0.0)).xyz
                 +1.0 * tex2D(ReShade::BackBuffer, texcoord + float2( d.x,  d.y)).xyz;

    float3 f_y = -1.0 * tex2D(ReShade::BackBuffer, texcoord + float2(-d.x, -d.y)).xyz
                 -2.0 * tex2D(ReShade::BackBuffer, texcoord + float2( 0.0, -d.y)).xyz
                 -1.0 * tex2D(ReShade::BackBuffer, texcoord + float2( d.x, -d.y)).xyz
                 +1.0 * tex2D(ReShade::BackBuffer, texcoord + float2(-d.x,  d.y)).xyz
                 +2.0 * tex2D(ReShade::BackBuffer, texcoord + float2( 0.0,  d.y)).xyz
                 +1.0 * tex2D(ReShade::BackBuffer, texcoord + float2( d.x,  d.y)).xyz;

    //Structure tensor components:
    float E = dot(f_x, f_x), F = dot(f_x, f_y), G = dot(f_y, f_y);
    
    //Eigenvalues
    float root = sqrt((E-G)*(E-G) + 4.0*F*F);
    float lambda1 = 0.5*(E+G + root);
    float lambda2 = 0.5*(E+G - root);

    //Eigenvector corresponding to larger eigenvalue
    float2 v = float2(lambda1 - E, -F);

    //Rotation matrix construction
    float arg = atan2(v.y, v.x);
    float s = sin(arg), c = cos(arg);
    float2x2 rot = float2x2(c, -s, s, c);

    //Scale transformation construction
    float ANI = (lambda1-lambda2)/(lambda1+lambda2);
    float2 scaling = float2(1.0+ANI, 1.0/(1.0+ANI));

    //Back to usual generalized filtering:
    const float PI  = 3.1415926535;
    const float GOLDEN_RATIO = 1.6180339887;

    float3 color = float3(0,0,0);

    //Number equivalent to standard technique
    const int samples = 4*(iKuwahara_Radius+1)*(iKuwahara_Radius+1);
    const float radius = float(iKuwahara_Radius);

    float3 means[8] = {float3(0,0,0), float3(0,0,0), float3(0,0,0), float3(0,0,0),
                       float3(0,0,0), float3(0,0,0), float3(0,0,0), float3(0,0,0)};

    float3 sigma[8] = {float3(0,0,0), float3(0,0,0), float3(0,0,0), float3(0,0,0),
                       float3(0,0,0), float3(0,0,0), float3(0,0,0), float3(0,0,0)};

    float norms[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    //Sample from disk, along a Fibbonaci spiral:
    const float2 scale = radius * BUFFER_PIXEL_SIZE;
    const float dt = 1.0/float(samples);
    const float A = 2.0*PI*float(samples)/GOLDEN_RATIO;

    float t = 0.0;

    for (int i=0; i<samples; i++) {
        //Generate point on the spiral:
        float ang = A*t;
        float2 offset = sqrt(t) * float2(sin(ang), cos(ang));

        //Take transformation from structure tensor into account:
        offset = scale * mul(rot, scaling*offset);

        //Retrieve index of the quadrant (angle is now also rotated)
        float f_ID = fmod(4.0*(ang-arg)/PI, 8.0);
        int i_ID = int(f_ID);

        //Gather values
        float3 val = tex2D(ReShade::BackBuffer, texcoord + offset).rgb;
                
        //For now this contains sums
        means[i_ID] += val;
        //And sums of squares
        sigma[i_ID] += val*val;
        //We now also need normalizations
        norms[i_ID] += 1.0;

        //Increment along the spiral
        t += dt;
    }

    //Large number
    float min_sigma = 1000.0;
    
    for (int k=0; k<8; k++) {
        //Now it's actually means
        means[k] /= norms[k];
        //And standard deviations
        sigma[k] = abs(sigma[k]/norms[k] - means[k]*means[k]);
        
        //Sum of standard deviations for each color channel
        float sigma2 = sigma[k].x + sigma[k].y + sigma[k].z;
        
        //Use mean from region with lowest deviation
        if (sigma2 < min_sigma) {
            min_sigma = sigma2;
            color = means[k];
        }
    }

    //=====Output:
	res.xyz = color;
	res.w = 1.0;
}

technique Kuwahara
{
	pass Kuwahara_Apply
	{
		VertexShader = PostProcessVS;
		PixelShader = PS_Kuwahara_Basic;
	}
}

technique KuwaharaGeneralized
{
	pass Kuwahara_Apply
	{
		VertexShader = PostProcessVS;
		PixelShader = PS_Kuwahara_Generalized;
	}
}

technique KuwaharaAnisotropic
{
    pass Kuwahara_Apply
    {
        VertexShader = PostProcessVS;
        PixelShader = PS_Kuwahara_Anisotropic;
    }
}