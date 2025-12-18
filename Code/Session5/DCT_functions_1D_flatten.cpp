#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using std::cout;
using std::endl;
using std::cos;
using std::string;

struct Sequence {
    int nbre_zero;  // number of zeros encountered
    std::vector<float> val;
};

const int width = 256;
const int height = 256;
const int vector_size = 256;
const int block = 8;

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

float *DCT_basis(size_t width_dim, size_t height_dim){
    float *DCT_pointer = new float[height_dim*width_dim];
    for (int k=0; k<height_dim; k++){
        for (int i =0; i<width_dim; i++){
            DCT_pointer[i+k*width_dim] = cos((k*M_PI/width_dim)*(i+0.5f));
        }
    }
    return DCT_pointer;
}

float *power(float *signal,size_t width_dim, size_t height_dim){
    float *power = new float[height_dim];
    for (int k=0; k<height_dim; k++){
        power[k]=0.0f;
        for (int i =0; i<width_dim; i++){
            power[k] += signal[i+k*width_dim]*signal[i+k*width_dim];
        }
    }
    return power;
}

float *normalization(float*signal, float* power,size_t width_dim, size_t height_dim){
    float *sign_norm = new float[height_dim*width_dim];
    for (int k=0; k<height_dim; k++){
        for (int i =0; i<width_dim; i++){
            sign_norm[i+k*width_dim] = signal[i+k*width_dim]/std::sqrt(power[k]);
        }
    }
    return sign_norm;
}

float *DCT_transform(float* image_vec,float *DCT_norm, size_t width_dim, size_t height_dim){
    float* DCT_coeff = new float [height_dim];
    for (int k=0; k<height_dim;k++){
        float sum = 0.0f;
        for (int i=0; i<width_dim;i++){
            sum +=DCT_norm[i+k*width_dim]*image_vec[i];
        }
        DCT_coeff[k] = sum;
    }
    return DCT_coeff;
}

float *Transpose_matrix(float* image, size_t width_dim, size_t height_dim){
    float* trans = new float[height_dim*width_dim];
    for (int k=0; k<height_dim; k++){
        for (int i=0; i<width_dim; i++){
            trans[k+i*width_dim]=image[i+k*width_dim];
        }
    }
    return trans;
}

float *Matrix_multiplicator(float * image_1, float * image_2, size_t width_dim, size_t height_dim){
    float *matrix_mu = new float[height_dim*width_dim];
    for (int w=0; w<height_dim; w++){
        for (int k=0; k<height_dim; k++){
            float sum = 0.0f;
            for (int i=0; i<width_dim; i++){
                sum += image_1[i+w*width_dim]*image_2[k+i*width_dim]; 
            }
            matrix_mu[k+w*width_dim] = sum;
        }
    }
    return matrix_mu;
}

float *Mat_vec_mul(float* image_vec, float*matrix, size_t width_dim, size_t height_dim){
    float* mat_vec_pointer = new float [height_dim*width_dim];
    for (int k=0; k<height_dim;k++){
        float sum = 0.0f;
        for (int i=0; i<width_dim;i++){
            sum +=matrix[i+k*width_dim]*image_vec[i];
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

float * JPEG_quantization(float * quant_table, float * image, size_t width_dim, size_t height_dim ){
    float* quantized_vector = new float[height_dim*width_dim];
    for (int j =0; j < height_dim; j++){
        for(int i = 0; i < width_dim; i++){
            
            quantized_vector[i+width_dim*j]= std::round(image[i+width_dim*j]/quant_table[i+width_dim*j]);
        }
    }
    return quantized_vector;

}

float * IJPEG(float * quant_table, float * image, size_t width_dim, size_t height_dim ){
    float* reconstructed_vector = new float[height_dim*width_dim];
    for (int j =0; j < height_dim; j++){
        for(int i = 0; i < width_dim; i++){
            
            reconstructed_vector[i+width_dim*j]= (image[i+width_dim*j]*quant_table[i+width_dim*j]);
        }
    }
    return reconstructed_vector;

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

float * Transform(float * Basis,float * image, size_t width_dim, size_t height_dim){
    float * DCT_matrix_1D = Matrix_multiplicator(Basis, image, width_dim, height_dim);
    float * trans_1D_DCT = Transpose_matrix(DCT_matrix_1D, width_dim, height_dim);
    float *DCT_matrix_2D= Matrix_multiplicator(Basis, trans_1D_DCT, width_dim, height_dim);
    float * trans_2D_DCT = Transpose_matrix(DCT_matrix_2D, width_dim, height_dim);

    delete [] DCT_matrix_1D;
    delete [] trans_1D_DCT;
    delete [] DCT_matrix_2D;
    return trans_2D_DCT;
}

float * IDCT_2D(float * image,float* Basis, size_t width_dim, size_t height_dim){

    float * IDCT_image = Matrix_multiplicator(image,Basis, width_dim, height_dim);
    float * trans_basis = Transpose_matrix(Basis, width_dim, height_dim);
    float * IDCT_image_2D= Matrix_multiplicator(trans_basis,IDCT_image, width_dim, height_dim);

    delete [] IDCT_image;
    delete [] trans_basis;
    return IDCT_image_2D;

}

void change_int_float(const int data[8][8],size_t height_dim,size_t width_dim, float * data_float){
    for (int j=0;j<height_dim;j++){
        for (int i=0;i<width_dim;i++){
            data_float[i+j*width_dim] = (float)data[j][i];
        }
    }
}

void approximate(float * image, float * quant_table, std::string file_2DCT, std::string file_2DCT_quant, std::string file_2DCT_IQ,std::string file_IDCT_quant, size_t height_dim, size_t width_dim, size_t block_size){
    int size = height_dim*width_dim;
    int block_elem = block_size*block_size;
    float * image_redist_block = new float [size];
    float * image_block = new float [block_elem];
    float * image_2DCT = new float [size];
    float * IQ_image_2DCT = new float[size];
    float * quant_image_2DCT = new float [size];
    float * quant_image_IDCT = new float [size];
    float * test_IDCT = new float [size];
    float * original_2DCT_image = new float [size];
    float * original_2DCT_quant = new float [size];
    float * original_2DCT_IQ = new float [size];
    float * original_IDCT_quant = new float [size];
    float * original = new float [size];


    int block_y = height/block_size;
    int block_x = width/block_size;


    if (size % block_elem == 0){
        int iteration = size/block_elem;
        float * DCT_basis_block = DCT_basis(block_size, block_size);
        float * power_block = power(DCT_basis_block,block_size,block_size);
        float * normal_DCT_base = normalization(DCT_basis_block, power_block,block_size, block_size);
        for (int k=0; k < block_y; k++){
            for (int j=0; j < block_x; j++){
                for (int w =0; w<block_size;w++){
                    for (int i=0 ; i<block_size ; i++ ){
                            image_redist_block[i+w*block_size+j*block_elem+(k*block_elem*block_x)] = image[(k*block_size*width_dim)+j*block_size+i+w*width_dim]; // Re build of vector with each group of 8 pixels aligned one after another instead that making a jump of 256 in the memory
                    }
                }
            }
        }
        for (int j=0; j<iteration;j++){
            for (int i=0 ; i<block_elem ; i++ ){
                image_block[i] = image_redist_block[i+j*block_elem];
                float *block_2DCT = Transform(normal_DCT_base,image_block,block_size,block_size);
                image_2DCT[i+j*block_elem] = block_2DCT [i];
                float *quant_2DCT = JPEG_quantization(quant_table,block_2DCT,block_size,block_size);
                quant_image_2DCT[i+j*block_elem] = quant_2DCT [i];
                float *IQ_2DCT = IJPEG(quant_table,quant_2DCT,block_size,block_size);
                IQ_image_2DCT[i+j*block_elem] = IQ_2DCT [i];
                float *block_IDCT = IDCT_2D(IQ_2DCT,normal_DCT_base,block_size,block_size);
                quant_image_IDCT[i+j*block_elem] = block_IDCT [i];

                float *test = IDCT_2D(block_2DCT,normal_DCT_base,block_size,block_size);
                test_IDCT[i+j*block_elem] = test [i];

            }
        }

        for (int k=0; k < block_y; k++){
            for (int j=0; j < block_x; j++){
                for (int w =0; w<block_size;w++){
                    for (int i=0 ; i<block_size ; i++ ){
                            original_2DCT_image[(k*block_size*width_dim)+j*block_size+i+w*width_dim] = image_2DCT[i+w*block_size+j*block_elem+(k*block_elem*block_x)];          // Re shape original matrice
                            original_2DCT_quant[(k*block_size*width_dim)+j*block_size+i+w*width_dim] = quant_image_2DCT[i+w*block_size+j*block_elem+(k*block_elem*block_x)];    // Re shape original matrice
                            original_2DCT_IQ[(k*block_size*width_dim)+j*block_size+i+w*width_dim] = IQ_image_2DCT[i+w*block_size+j*block_elem+(k*block_elem*block_x)];          // Re shape original matrice
                            original_IDCT_quant[(k*block_size*width_dim)+j*block_size+i+w*width_dim] = quant_image_IDCT[i+w*block_size+j*block_elem+(k*block_elem*block_x)];    // Re shape original matrice
                            original[(k*block_size*width_dim)+j*block_size+i+w*width_dim] = test_IDCT[i+w*block_size+j*block_elem+(k*block_elem*block_x)];                      // Re shape original matrice
                    }
                }
            }
        }
        
        Raw_save(file_2DCT,original_2DCT_image,size);
        Raw_save(file_2DCT_quant,original_2DCT_quant,size);
        Raw_save(file_2DCT_IQ,original_2DCT_IQ,size);
        Raw_save(file_IDCT_quant,original_IDCT_quant,size);
        Raw_save("Original_IDCT.raw",original,size);

        delete [] DCT_basis_block;
        delete [] power_block;
        delete [] normal_DCT_base;
    }
    else {
        cout << "ERROR - NOT DIVISIBLE IN BLOCK" <<endl;
    }

    delete [] image_redist_block;
    delete [] image_block;
    delete [] image_2DCT;
    delete [] IQ_image_2DCT;
    delete [] quant_image_2DCT;
    delete [] quant_image_IDCT;
    delete [] test_IDCT;
    delete [] original_2DCT_image;
    delete [] original_2DCT_quant;
    delete [] original_2DCT_IQ;
    delete [] original_IDCT_quant;
    delete [] original;
}

float* encoder(float * image, float * quant_table, size_t height_dim, size_t width_dim, size_t block_size){

    int size = height_dim*width_dim;
    int block_x,block_y;

    float* quantized_image = new float [size];
    float* image_block = new float [size];

    block_y = height_dim/block_size;
    block_x = width_dim/block_size;

    float * DCT_basis_block = DCT_basis(block_size, block_size);
    float * power_block = power(DCT_basis_block,block_size,block_size);
    float * normal_DCT_base = normalization(DCT_basis_block, power_block,block_size, block_size);

    for (int by=0; by < block_y; by++){
        for (int bx=0 ; bx < block_x; bx++){
            for (int j=0; j < block_size; j++){
                for (int i = 0; i < block_size; i++){
                    image_block[i+block_size*j] = image[(by*8 + j)*width_dim + (bx*8 + i)];
                }
            }

            float* DCT_block = Transform(normal_DCT_base, image_block, 8, 8);

            // --- Quantization ---
            float* Quantized_block = JPEG_quantization(quant_table, DCT_block, 8, 8);

            for (int j=0; j < block_size; j++){
                for (int i = 0; i < block_size; i++){
                    quantized_image[(by*8 + j)*width_dim + (bx*8 + i)] = Quantized_block[i+block_size*j];
                }
            }

            delete [] DCT_block;
            delete [] Quantized_block;
        }
    }

    delete [] image_block;
    delete [] DCT_basis_block;
    delete [] power_block;
    delete [] normal_DCT_base;
    //return blocks;
    return quantized_image;
} 

float* decoder(float * quantized_image, float * quant_table, size_t height_dim, size_t width_dim, size_t block_size){

    int size = height_dim*width_dim;
    int block_elem = block_size*block_size;
    int block_x,block_y;

    float* image = new float [size];
    float* IDCT_block = new float [size];
    float* extract_quantized = new float [block_elem];

    block_y = height_dim/block_size;
    block_x = width_dim/block_size;

    float * DCT_basis_block = DCT_basis(block_size, block_size);
    float * power_block = power(DCT_basis_block,block_size,block_size);
    float * normal_DCT_base = normalization(DCT_basis_block, power_block,block_size, block_size);

    int b = 0;
    float dequant_block[block_elem];

    for (int by=0; by < block_y; by++){
        for (int bx=0 ; bx < block_x; bx++){
            for (int j=0; j < block_size; j++){
                for (int i = 0; i < block_size; i++){
                    extract_quantized[i+block_size*j] = quantized_image[(by*8 + j)*width_dim + (bx*8 + i)];
                }
            }
            
            for (int k=0; k < block_elem; k++){
                dequant_block[k] = extract_quantized[k]*quant_table[k];
            }

            IDCT_block = IDCT_2D(dequant_block,normal_DCT_base,block_size,block_size);

            for (int j=0; j < block_size; j++){
                for (int i = 0; i < block_size; i++){
                    image[(by*8 + j)*width_dim + (bx*8 + i)] = IDCT_block[i+block_size*j];
                }
            }

            delete [] IDCT_block;
        }
    }
    delete [] DCT_basis_block;
    delete [] power_block;
    delete [] normal_DCT_base;
    return image;
} 

float* DC_encoder(float * image, float * quant_table, size_t height_dim, size_t width_dim, size_t block_size){

    int size = height_dim*width_dim;
    int block_elem = block_size*block_size;
    int block_x,block_y;
    int DC_component = size/block_elem;

    float* DC_block = new float [DC_component];
    float* image_block = new float [size];


    block_y = height_dim/block_size;
    block_x = width_dim/block_size;

    float * DCT_basis_block = DCT_basis(block_size, block_size);
    float * power_block = power(DCT_basis_block,block_size,block_size);
    float * normal_DCT_base = normalization(DCT_basis_block, power_block,block_size, block_size);

    int k=0;

    for (int by=0; by < block_y; by++){
        for (int bx=0 ; bx < block_x; bx++){
            for (int j=0; j < block_size; j++){
                for (int i = 0; i < block_size; i++){
                    image_block[i+block_size*j] = image[(by*8 + j)*width_dim + (bx*8 + i)];
                }
            }

            float* DCT_block = Transform(normal_DCT_base, image_block, 8, 8);

            // --- Quantization ---
            float* Quantized_block = JPEG_quantization(quant_table, DCT_block, 8, 8);

            DC_block[k] = Quantized_block[0];
            k++;

            delete [] DCT_block;
            delete [] Quantized_block;
        }
    }
    delete [] image_block;
    delete [] DCT_basis_block;
    delete [] power_block;
    delete [] normal_DCT_base;
    return DC_block;
}

void difference_text(const string& filename, float* DC_coef, size_t DC_size){

    std::ofstream txt(filename);

    float prev = 0;
    for (size_t i = 0; i < DC_size; i++) {
        float dc = DC_coef[i];
        float delta = dc-prev;
        txt << delta << "\n";
        prev = dc;
    }
    txt.close();
}

float* reconstruct_from_txt(const string& filename,size_t vec_size){
    float* reconstruct_DC = new float [vec_size];
    std::ifstream in(filename);

    int diff;
    int prev = 0;
    int k = 0;
    bool first = true;

    while (in >> diff) {
        int dc;
        if (first) {
            dc = diff;
            first = false;
        } else {
            dc = prev + diff;
        }
        reconstruct_DC[k] = dc;
        prev = dc;
        k++;
    }
    return reconstruct_DC;
}

std::ostream& operator<<(std::ostream& os, const Sequence& seq)
{
    os << "{ zeros = " << seq.nbre_zero << ", values = [";
    for (size_t i = 0; i < seq.val.size(); i++)
    {
        os << seq.val[i];
        if (i + 1 < seq.val.size())
            os << ", ";
    }
    os << "] }";
    return os;
}

std::vector<Sequence> AC_sequence(float* quantized_block, size_t vec_size){
    
    std::vector<Sequence> Full_seq;
    int i = 1; // begin at index 1 - first AC element

    while ( i < vec_size){
        int nbre_zero = 0;
        while ( i < vec_size && quantized_block[i]==0){
            nbre_zero ++;
            i ++;
        }

        Sequence seq;
        seq.nbre_zero = nbre_zero;

        while ( i < vec_size && quantized_block[i]!=0){
            seq.val.push_back(quantized_block[i]);
            i ++; 
        } 
        Full_seq.push_back(seq);
    }
    return Full_seq;
}

float* AC_sequence_decode(float dc, std::vector<Sequence> AC_encoded, size_t vec_size){

    size_t list_size = AC_encoded.size();
    float* AC_reconstr = new float [vec_size];
    int k = 0;
    AC_reconstr[k] = dc;
    k++;

    for (int i=0 ; i < list_size ; i++){
        int zeros = AC_encoded[i].nbre_zero;
        for (int j = 0; j < zeros; j++){
            AC_reconstr[k] = 0;
            k++;
        }
        int AC_size = AC_encoded[i].val.size();
        for (int j = 0; j < AC_size; j++){
            AC_reconstr[k] = AC_encoded[i].val[j];
            k++;
        }
    }
    return AC_reconstr;
}

std::vector <std::vector<Sequence>> encoder_RLE(float * image, float * quant_table, float * dc_coef, size_t height_dim, size_t width_dim, size_t block_size){

    int size = height_dim*width_dim;
    int block_elem = block_size*block_size;
    int block_x,block_y;
    int visu = 0;

    std::vector <std::vector<Sequence>> quantized_image;
    float* image_block = new float [size];

    block_y = height_dim/block_size;
    block_x = width_dim/block_size;

    float * DCT_basis_block = DCT_basis(block_size, block_size);
    float * power_block = power(DCT_basis_block,block_size,block_size);
    float * normal_DCT_base = normalization(DCT_basis_block, power_block,block_size, block_size);

    for (int by=0; by < block_y; by++){
        for (int bx=0 ; bx < block_x; bx++){
            for (int j=0; j < block_size; j++){
                for (int i = 0; i < block_size; i++){
                    image_block[i+block_size*j] = image[(by*8 + j)*width_dim + (bx*8 + i)];
                }
            }

            float* DCT_block = Transform(normal_DCT_base, image_block, 8, 8);

            // *** Quantization ***
            float* Quantized_block = JPEG_quantization(quant_table, DCT_block, 8, 8);

            dc_coef[by*block_x+bx] = Quantized_block[0];

            std::vector<Sequence> AC_encod = AC_sequence(Quantized_block,block_elem);

            if (visu == 4){
                size_t list_encod = AC_encod.size();
                for (int i =0 ; i < list_encod ; i++){
                    cout << "RLE sequence for the "<< visu << "th block of 8x8" << AC_encod[i] << endl;
                }
            }
            visu ++;

            quantized_image.push_back(AC_encod);

            delete [] DCT_block;
            delete [] Quantized_block;
        }
    }

    delete [] image_block;
    delete [] DCT_basis_block;
    delete [] power_block;
    delete [] normal_DCT_base;
    //return blocks;
    return quantized_image;
} 

float* decoder_RLE(float * quant_table, std::vector <std::vector<Sequence>> AC, float* dc_coef, size_t height_dim, size_t width_dim, size_t block_size){

    int size = height_dim*width_dim;
    int block_elem = block_size*block_size;
    int block_x,block_y;
    bool visu = true;

    std::vector <std::vector<Sequence>> quantized_image;
    float* image_block = new float [size];
    float* image_return = new float [size];
    float* IDCT_block = new float [size];
    float* extract_quantized = new float [block_elem];
    float* dequant_block = new float [block_elem];
    float* reconstruct_DCT = new float [size];

    block_y = height_dim/block_size;
    block_x = width_dim/block_size;

    float * DCT_basis_block = DCT_basis(block_size, block_size);
    float * power_block = power(DCT_basis_block,block_size,block_size);
    float * normal_DCT_base = normalization(DCT_basis_block, power_block,block_size, block_size);

    size_t list_size = AC.size();
    cout << "Size of AC RLE :" << list_size <<endl;
    for (int i = 0; i < list_size ; i++){
        float* reconstruct_block = AC_sequence_decode(dc_coef[i], AC[i],block_elem);
        for (int k=0; k < block_elem; k++){
            reconstruct_DCT[i*block_elem+k] = reconstruct_block[k];
        }
    }

    for (int by=0; by < block_y; by++){
        for (int bx=0 ; bx < block_x; bx++){
            for (int j=0; j < block_size; j++){
                for (int i = 0; i < block_size; i++){
                    extract_quantized[i+block_size*j] = reconstruct_DCT[(by*block_x+bx)*block_elem+i+block_size*j];
                }
            }
        

            for (int k=0; k < block_elem; k++){
                dequant_block[k] = extract_quantized[k]*quant_table[k];
            }

            IDCT_block = IDCT_2D(dequant_block,normal_DCT_base,block_size,block_size);
    
            for (int j=0; j < block_size; j++){
                for (int i = 0; i < block_size; i++){
                    image_return[(by*8 + j)*width_dim + (bx*8 + i)] = IDCT_block[i+block_size*j];
                }
            }
            delete [] IDCT_block;
            
        }

    }

    delete [] image_block;
    delete [] DCT_basis_block;
    delete [] power_block;
    delete [] normal_DCT_base;
    //return blocks;
    return image_return;
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
    float*DCT_basis_pointer = DCT_basis(width,height);
    Raw_save("DCT_basis.raw",DCT_basis_pointer,vector_size*vector_size);
    x = 125;
    y = 125;
    cout << "Check DCT basis content at row " << y << " and column " << x << " = " << DCT_basis_pointer[x+width*y] << endl;
    float * power_DCT_basis = power(DCT_basis_pointer,width,height);
    float * DCT_basis_normalization = normalization(DCT_basis_pointer,power_DCT_basis,width,height);
    //float* DCT_coeff = DCT_transform(parrot_loading,width);
    //cout << "Check DCT coeff for first vector of parrot at position " << y << " is " << DCT_coeff[y] << endl;
    float* DCT_basis_trans = Transpose_matrix(DCT_basis_normalization,width,height);
    x=5;
    y=10;
    cout << "Normed of DCT at row " << y << " and column " << x << " = " << DCT_basis_normalization[x+width*y]<< endl;
    x=10;
    y=5;
    cout << "Transposed of DCT at row " << y << " and column " << x << " = " << DCT_basis_trans[x+width*y]<< endl;
    float * ortho_matrix = Matrix_multiplicator(DCT_basis_normalization,DCT_basis_trans, width,height);
    x = 20;
    y = 20;
    cout << "Ortho check for row " << y << " and column " << x << " is " << ortho_matrix[x+width*y] << endl;

    // Parrot DCT Transform
    float* DCT_parrot_1D = DCT_transform(parrot_loading,DCT_basis_normalization, width, height);
    Raw_save("DCT_parrot.raw",DCT_parrot_1D,vector_size);
    // Cosine DCT transform
    float* DCT_cosine = DCT_transform(&cosine_raw[0][0],DCT_basis_normalization, width, height);
    Raw_save("DCT_cosine.raw",DCT_cosine,vector_size);
    // Noise DCT transofrm
    float* noise_pointer = Load_raw("Noise_8_std.raw",height,width);
    float* DCT_noise = DCT_transform(noise_pointer,DCT_basis_normalization,width, height);
    Raw_save("DCT_noise.raw",DCT_noise,vector_size);

    // Quantization
    int threshold_val = 100;
    float* Quantized_parrot_1D = Quantization(threshold_val,DCT_parrot_1D,width);
    int parrot_count_zero = count_zero(Quantized_parrot_1D, width);
    cout <<"Nbre of zeros in quantized parrot : " << parrot_count_zero << endl;
    Raw_save("Quantized_parrot.raw",Quantized_parrot_1D,vector_size);
    float* IDCT_parrot = Mat_vec_mul(Quantized_parrot_1D,DCT_basis_trans,width, height);
    float parrot_psnr = psnr(IDCT_parrot,parrot_loading,width);
    cout <<"PSNR in quantized parrot with threshold of "<<threshold_val<<" : " << parrot_psnr << " dB" << endl;

    // ************************ SESSION 4 *************************

    // 2DC TRANSFORM
    float* DCT_parrot_2D = Transform(DCT_basis_normalization,parrot_loading, width, height);
    Raw_save("2DCT_parrot.raw",DCT_parrot_2D,height*width);
    int threshold_val_2D = 50;
    float* Quantized_parrot_2D = Quantization(threshold_val_2D,DCT_parrot_2D,height*width);
    int parrot_count_zero_2D = count_zero(Quantized_parrot_2D, height*width);
    cout <<"Nbre of zeros in quantized parrot : " << parrot_count_zero_2D << endl;
    Raw_save("2DCT_Quantized_parrot.raw",Quantized_parrot_2D,vector_size);
    float * IDCT_parrot_2D = IDCT_2D(Quantized_parrot_2D,DCT_basis_normalization, width, height);
    float parrot_psnr_2D = psnr(IDCT_parrot_2D,parrot_loading,height*width);
    Raw_save("IDCT_Quantized_parrot_2D.raw",IDCT_parrot_2D,height*width);
    cout <<"PSNR in 2D Quantized Parrot with threshold of "<<threshold_val_2D<<" : " << parrot_psnr_2D << " dB" << endl;

    // JPEG QUANTIZATION - 50% Quality

    const int JPEG_quant[8][8] = {
    {16, 11, 10, 16, 24, 40,  51,  61},
    {12, 12, 14, 19, 26, 58,  60,  55},
    {14, 13, 16, 24, 40, 57,  69,  56},
    {14, 17, 22, 29, 51, 87,  80,  62},
    {18, 22, 37, 56, 68, 109, 103, 77},
    {24, 35, 55, 64, 81, 104, 113, 92},
    {49, 64, 78, 87, 103, 121, 120, 101},
    {72, 92, 95, 98, 112, 100, 103, 99}
    };

    float Q_float[64];
    change_int_float(JPEG_quant, 8, 8, Q_float);
    Raw_save("JPEG_50.raw",Q_float,64);
    
    // 8x8 Block

    approximate(parrot_loading,Q_float,"2DCT_parrot_block.raw","2DCT_parrot_block_quant.raw","2DCT_parrot_block_IQ.raw","IDCT_parrot_block_quant.raw",height,width,block);


    // ************************ SESSION 5 *************************
    float* encoded_data = encoder(parrot_loading,Q_float,height,width,block);
    Raw_save("Encoded_data.raw",encoded_data,height*width);
    float* decoder_data = decoder(encoded_data,Q_float,height,width,block);
    Raw_save("Decoded_data.raw",decoder_data,height*width);
    float* DC_quant= DC_encoder(parrot_loading,Q_float,height,width,block);
    Raw_save("DC_quant.raw",DC_quant,1024);

    difference_text("DC_coeff_difference.txt",DC_quant,1024);           
    // Difference very close to 0 except when we jump to another line (each 32 elements) and with first element
    // Idea - use a tracking in zig-zag to increase spatial correlation and reduce difference for minimum encoding
    // Second idea - average the elements over the mean of all DC to reduce the values.

    float* reconstr_DC_coef = reconstruct_from_txt("DC_coeff_difference.txt",1024);
    Raw_save("reconstr_DC_coef.raw",reconstr_DC_coef,1024);

    float* dc_coef = new float [height*width/(block*block)];
    std::vector<std::vector<Sequence>> Encoded_RLE_data = encoder_RLE(parrot_loading,Q_float,dc_coef,height,width,block);
    // Improve the encoding with euclidian distance sorting !
    float* image_decode_after_RLE = decoder_RLE(Q_float,Encoded_RLE_data,dc_coef,height,width,block);
    Raw_save("reconstr_after_RLE.raw",image_decode_after_RLE,height*width);

    // ************************ SESSION 6 *************************






    // ************************ MEMORY FREE *************************
    delete [] cosine_raw;
    delete [] parrot_loading;
    delete [] multiplied_parrot;
    delete [] DCT_basis_pointer;
    //delete [] DCT_coeff;
    delete [] DCT_basis_trans;
    delete [] ortho_matrix;
    delete [] noise_pointer;
    delete [] DCT_parrot_1D;
    delete [] DCT_cosine;
    delete [] DCT_noise;
    delete [] Quantized_parrot_1D;
    delete [] IDCT_parrot;
    delete [] DCT_parrot_2D;
    delete [] Quantized_parrot_2D;
    delete [] IDCT_parrot_2D;
    delete [] encoded_data;
    delete [] decoder_data;
    delete [] dc_coef;
    delete [] image_decode_after_RLE;

    return EXIT_SUCCESS;
}


