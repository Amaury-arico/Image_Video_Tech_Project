#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using std::cout;
using std::endl;
using std::cos;
using std::string;

const int width = 256;
const int height = 256;
const int vector_size = 256;

float (*cosine_image(int h))[width]{
    // dynamic allocation
    float (*image)[width] = new float[h][width];
    // loop for image content
    for (int y=0; y<h; y++){
        for (int x=0; x<width; x++){
            image[y][x] = 0.5f+0.5f*cos(x*M_PI/32.0f)*cos(y*M_PI/64.0f);
        }
    }
    // return pointer
    return image;
}

void Raw_save(const string& filename, const float* data, size_t size){
    std::ofstream out(filename, std::ios::binary);

    if(!out){
        std::cerr << "Error ! coudn't open the file " << filename << endl;
        return;
    }
    //float* temp = &data[0];
    char *data_point = (char*) data;
    out.write(data_point, size*sizeof(float));
    out.close();
}

float *Load_raw(const string& filename, const int h, const int w){
    std::ifstream in(filename, std::ios::binary);

    float* data = new float[h*w];
    if(!in){
        std::cerr << "Cannot open the file " << filename << endl;
        return nullptr;
    }
    char* byte_pointer = (char*) data;
    in.read(byte_pointer, h*w*sizeof(float));
    return data;
}

float *multiplicator (float* image_1, float* image_2, size_t size){

    float* multiplied_image = new float [size];
    for (int i=0; i<size; i++){
        multiplied_image[i] = image_1[i]*image_2[i];
    }
    return multiplied_image;
}

float mse(float* image_1, float* image_2, size_t size){
    float mse_value = 0.0f;
    float sum = 0.0f;
    float diff = 0.0f;
    for (int i=0; i<size; i++ ){
        diff = image_1[i] - image_2[i];
        sum += diff * diff;
    }
    mse_value = sum/size;
    return mse_value;
}

float max_val(float* image, size_t size){
    float max = image[0];
    for (int i=0; i<size; i++){
        if (max < image[i]){
            max = image[i];
        }
    }
    return max;
}

float psnr(float* image_1, float* image_2, size_t size){
    float psnr_val = 0.0f;
    float mse_val = mse(image_1, image_2, size);
    float max = max_val(image_1,size);

    psnr_val = 10.0f*std::log10(max*max/mse_val);
    return psnr_val;
}

float (*DCT_basis())[vector_size]{
    float (*DCT_pointer)[vector_size] = new float[vector_size][vector_size];
    for (int k=0; k<vector_size; k++){
        for (int i =0; i<vector_size; i++){
            DCT_pointer[k][i] = cos((k*M_PI/vector_size)*(i+0.5f));
        }
    }
    return DCT_pointer;
}

float *power(float (*signal)[vector_size]){
    float *power = new float[vector_size];
    for (int k=0; k<vector_size; k++){
        power[k]=0.0f;
        for (int i =0; i<vector_size; i++){
            power[k] += signal[k][i]*signal[k][i];
        }
    }
    return power;
}

float (*normalization(float(*signal)[vector_size], float* power))[vector_size]{
    float (*sign_norm)[vector_size] = new float[vector_size][vector_size];
    for (int k=0; k<vector_size; k++){
        for (int i =0; i<vector_size; i++){
            sign_norm[k][i] = signal[k][i]/std::sqrt(power[k]);
        }
    }
    return sign_norm;
}

float *DCT_transform(float* image_vec,float(*DCT_norm)[vector_size], size_t size){
    float* DCT_coeff = new float [size];
    for (int k=0; k<vector_size;k++){
        float sum = 0.0f;
        for (int i=0; i<size;i++){
            sum +=DCT_norm[k][i]*image_vec[i];
        }
        DCT_coeff[k] = sum;
    }
    return DCT_coeff;
}

float (*Transpose_matrix(float(* image)[vector_size], const int h))[vector_size]{
    float(* trans)[vector_size] = new float[h][vector_size];
    for (int k=0; k<h; k++){
        for (int i=0; i<vector_size; i++){
            trans[k][i]=image[i][k];
        }
    }
    return trans;
}

float(*Matrix_multiplicator(float(* image_1)[vector_size], float(* image_2)[vector_size]))[vector_size]{
    float (*matrix_mu)[vector_size] = new float[vector_size][vector_size];
    for (int w=0; w<vector_size; w++){
        for (int k=0; k<vector_size; k++){
            float sum = 0.0f;
            for (int i=0; i<vector_size; i++){
                sum += image_1[w][i]*image_2[i][k]; 
            }
            matrix_mu[w][k] = sum;
        }
    }
    return matrix_mu;
}

float *Mat_vec_mul(float* image_vec, float(*matrix)[vector_size], size_t size){
    float* mat_vec_pointer = new float [size];
    for (int k=0; k<vector_size;k++){
        float sum = 0.0f;
        for (int i=0; i<size;i++){
            sum +=matrix[k][i]*image_vec[i];
        }
        mat_vec_pointer[k] = sum;
    }
    return mat_vec_pointer;
}

float* Quantization(float threshold,float* image,size_t size ){
    float* quantized_vector = new float[size];
    for(int i=0;i<size;i++){
        if(std::abs(image[i])>threshold){
            quantized_vector[i]=image[i];
        }
        else{
            quantized_vector[i]=0.0f;
        }
    }
    return quantized_vector;
}

int count_zero(float* image,size_t size ){
    int count = 0;
    for(int i=0;i<size;i++){
        if(std::abs(image[i])<1e-6f){
            count ++;
        }
    }
    return count;
}

int main(){

    // ************************ SESSION 2 *************************
    // bi-dimensional image pointer
    float (*cosine_raw)[width] = cosine_image(height);

    // printf
    int x = 1;
    int y = 1;

    cout << "Check image content at row " << y << "and column " << x << " = " << cosine_raw[y][x] << endl;
    Raw_save("output_grey_raw.raw",&cosine_raw[0][0],height*width);
    float *parrot_loading = Load_raw("parrot_256x256.raw",height,width);
    float *multiplied_parrot = multiplicator(parrot_loading,&cosine_raw[0][0],height*width);
    Raw_save("parrot_cosine.raw", multiplied_parrot, height*width);
    float mse_val = mse(multiplied_parrot, parrot_loading, height*width);
    float max = max_val(multiplied_parrot, height*width);
    float psnr_val = psnr(multiplied_parrot, parrot_loading, height*width);

    //cout << "value of cosined parrot [1] is "<< multiplied_parrot[1] << endl;
    //cout << "value of parrot [1] is "<< parrot_loading[1] << endl;
    //cout << "MSE of cosined parrot is "<< mse_val << endl;
    //cout << "Max of cosined parrot is "<< max << endl;
    //cout << "PSNR of cosined parrot is "<< psnr_val << endl;

    // ************************ SESSION 3 *************************
    float(*DCT_basis_pointer)[vector_size] = DCT_basis();
    Raw_save("DCT_basis.raw",&DCT_basis_pointer[0][0],vector_size*vector_size);
    x = 125;
    y = 125;
    cout << "Check DCT basis content at row " << y << " and column " << x << " = " << DCT_basis_pointer[y][x] << endl;
    float* power_DCT_basis = power(DCT_basis_pointer);
    float (*DCT_basis_normalization)[vector_size] = normalization(DCT_basis_pointer,power_DCT_basis);
    //float* DCT_coeff = DCT_transform(parrot_loading,width);
    //cout << "Check DCT coeff for first vector of parrot at position " << y << " is " << DCT_coeff[y] << endl;
    float(* DCT_basis_trans)[vector_size] = Transpose_matrix(DCT_basis_normalization,vector_size);
    x=5;
    y=10;
    cout << "Normed of DCT at row " << y << " and column " << x << " = " << DCT_basis_normalization[y][x]<< endl;
    x=10;
    y=5;
    cout << "Transposed of DCT at row " << y << " and column " << x << " = " << DCT_basis_trans[y][x]<< endl;
    float(* ortho_matrix)[vector_size] = Matrix_multiplicator(DCT_basis_normalization,DCT_basis_trans);
    x = 20;
    y = 20;
    cout << "Ortho check for row " << y << " and column " << x << " is " << ortho_matrix[y][x] << endl;

    // Parrot DCT Transform
    float* DCT_parrot = DCT_transform(parrot_loading,DCT_basis_normalization,width);
    Raw_save("DCT_parrot.raw",DCT_parrot,vector_size);
    // Cosine DCT transform
    float* DCT_cosine = DCT_transform(&cosine_raw[0][0],DCT_basis_normalization,width);
    Raw_save("DCT_cosine.raw",DCT_cosine,vector_size);
    // Noise DCT transofrm
    float* noise_pointer = Load_raw("Noise_8_std.raw",height,width);
    float* DCT_noise = DCT_transform(noise_pointer,DCT_basis_normalization,width);
    Raw_save("DCT_noise.raw",DCT_noise,vector_size);

    // Quantization
    int threshold_val = 100;
    float* Quantized_parrot = Quantization(threshold_val,DCT_parrot,width);
    int parrot_count_zero = count_zero(Quantized_parrot, width);
    cout <<"Nbre of zeros in quantized parrot : " << parrot_count_zero << endl;
    Raw_save("Quantized_parrot.raw",Quantized_parrot,vector_size);
    float* IDCT_parrot = Mat_vec_mul(Quantized_parrot,DCT_basis_trans,width);
    float parrot_psnr = psnr(IDCT_parrot,parrot_loading,width);
    cout <<"PSNR in quantized parrot with threshold of "<<threshold_val<<" : " << parrot_psnr << " dB" << endl;


    // ************************ MEMORY FREE *************************
    delete [] cosine_raw;
    delete [] parrot_loading;
    delete [] multiplied_parrot;
    delete [] DCT_basis_pointer;
    //delete [] DCT_coeff;
    delete [] DCT_basis_trans;
    delete [] ortho_matrix;
    delete [] noise_pointer;
    delete [] DCT_parrot;
    delete [] DCT_cosine;
    delete [] DCT_noise;
    delete [] Quantized_parrot;
    delete [] IDCT_parrot;

    return EXIT_SUCCESS;
}


