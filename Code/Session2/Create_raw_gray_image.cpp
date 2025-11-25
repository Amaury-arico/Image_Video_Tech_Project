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

float (*grey_image(int h))[width]{
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

int main(){
    // bi-dimensional image pointer
    float (*raw_image)[width] = grey_image(height);

    // printf
    int x = 1;
    int y = 1;

    cout << "Check image content at row " << y << "and column " << x << " = " << raw_image[y][x] << endl;
    Raw_save("output_grey_raw.raw",&raw_image[0][0],height*width);
    float *parrot_loading = Load_raw("parrot_256x256.raw",height,width);
    float *multiplied_parrot = multiplicator(parrot_loading,&raw_image[0][0],height*width);
    Raw_save("parrot_cosine.raw", multiplied_parrot, height*width);
    float mse_val = mse(multiplied_parrot, parrot_loading, height*width);
    float max = max_val(multiplied_parrot, height*width);
    float psnr_val = psnr(multiplied_parrot, parrot_loading, height*width);

    cout << "value of cosined parrot [1] is "<< multiplied_parrot[1] << endl;
    cout << "value of parrot [1] is "<< parrot_loading[1] << endl;
    cout << "MSE of cosined parrot is "<< mse_val << endl;
    cout << "Max of cosined parrot is "<< max << endl;
    cout << "PSNR of cosined parrot is "<< psnr_val << endl;


    delete [] multiplied_parrot;
    delete [] parrot_loading;
    delete [] raw_image;

    return EXIT_SUCCESS;
}


