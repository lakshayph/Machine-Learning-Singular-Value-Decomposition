function svd_power( training_file, m, iterations)
training_data  = importdata(training_file);
size_projection = str2double(m);
num_iterations = str2double(iterations);
storage_u = zeros(size_projection,size(training_data,1));
storage_v = zeros(size_projection,size(training_data,2));
s_matrix = zeros(size_projection,size_projection);

initial_matrix_u = training_data * training_data';
initial_matrix_v = training_data' * training_data;

    for j=1:size_projection
        b_u = ones(size(initial_matrix_u,2),1);
        b_v = ones(size(initial_matrix_v,2),1);
        u = initial_matrix_u;
        v = initial_matrix_v;
        for k =1:num_iterations
            a_u = u * b_u;
            c_u = norm(a_u);
            b_u = a_u/c_u;
            
            a_v = v * b_v;
            c_v = norm(a_v);
            b_v = a_v/c_v;
        end
        dummy_u = b_u;
        dummy_v = b_v;
        s_matrix(j,j) = sqrt((dummy_u' * u) * dummy_u);
        initial_matrix_u = initial_matrix_u - ((initial_matrix_u*dummy_u)*dummy_u');
        initial_matrix_v = initial_matrix_v - ((initial_matrix_v*dummy_v)*dummy_v');
        storage_u(j,:) = dummy_u(:,1);
        storage_v(j,:) = dummy_v(:,1);       
    end
    
    svd = (storage_u' * s_matrix) * storage_v;
        
    fprintf('Matrix U: \n')
    for i = 1:size(training_data,1)
         fprintf('Row %3d:',i);
        for j =1:size_projection
            fprintf('%8.4f', storage_u(j,i));
        end
        fprintf('\n');
    end
    
    fprintf('\nMatrix S: \n')
    for i = 1:size_projection
        fprintf('Row %3d:',i);
        for j =1:size_projection
            fprintf('%8.4f', s_matrix(i,j));
        end
        fprintf('\n');
    end
    
    fprintf('\nMatrix V: \n')
    for i = 1:size(training_data,2)
         fprintf('Row %3d:',i);
        for j =1:size_projection
            fprintf('%8.4f', storage_v(j,i));
        end
        fprintf('\n');
    end
    
    fprintf('\nReconstruction:\n')
    for i = 1:size(training_data,1)
         fprintf('Row %3d:',i);
        for j =1:size(training_data,2)
            fprintf('%8.4f', svd(i,j));
        end
        fprintf('\n');
    end    
end