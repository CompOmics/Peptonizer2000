use serde::{Serialize, Deserialize };

pub trait HttpClient {
    fn perform_post_request<T: Serialize>(&self, url: String, batch: &T) -> Result<String, String>;
}

#[cfg(target_arch = "wasm32")]
pub struct WasmHttpClient;

#[cfg(not(target_arch = "wasm32"))]
pub struct PyHttpClient;

#[cfg(target_arch = "wasm32")]
impl HttpClient for WasmHttpClient {

    fn perform_post_request<T: Serialize>(&self, url: String, payload: &T) -> Result<String, String> {
        use web_sys::{XmlHttpRequest};

        let payload_json = serde_json::to_string(payload).unwrap();

        // Create a new XMLHttpRequest object
        let xhr = XmlHttpRequest::new().map_err(|e| format!("Failed to create XMLHttpRequest: {:?}", e))?;
        
        // Open the request (synchronous mode by setting `async` to false)
        xhr.open_with_async("POST", &url, false)
            .map_err(|e| format!("Failed to open request: {:?}", e))?;
        
        // Set the request header for JSON
        xhr.set_request_header("Content-Type", "application/json")
            .map_err(|e| format!("Failed to set request header: {:?}", e))?;
        
        // Send the request with the body
        xhr.send_with_opt_str(Some(&payload_json))
            .map_err(|e| format!("Failed to send request: {:?}", e))?;
        
        let status = xhr.status().map_err(|_e| format!("Failed to extract status from response"))?;
        if status == 200 {
            let response = xhr.response_text()
            .expect("Expected json in response")
            .ok_or(format!("Failed to extract text from response"))?;
        
            return Ok(format!("{}", response));
        }

        Err(format!("Status code {}", status))
    }
}

#[cfg(not(target_arch = "wasm32"))]
impl HttpClient for PyHttpClient {

    fn perform_post_request<T: Serialize>(&self, url: String, payload: &T) -> Result<String, String> {
        use reqwest::Client;
        use tokio::runtime::Runtime;

        // Create a Tokio runtime for async execution
        let rt = Runtime::new().unwrap();

        // Execute the HTTP POST request within the runtime
        let result = rt.block_on(async {
            let client = Client::new();
            let response = client.post(&url)
                .json(payload)
                .send()
                .await?;

            // Get the response body as a string
            response.text().await
        });
        
        // Handle the result and convert to PyResult
        match result {
            Ok(body) => Ok(body),
            Err(e) => Err(format!("HTTP POST request failed: {}", e)),
        }
    }
}

#[cfg(target_arch = "wasm32")]
pub fn create_http_client() -> impl HttpClient {
    WasmHttpClient
}

#[cfg(not(target_arch = "wasm32"))]
pub fn create_http_client() -> impl HttpClient {
    PyHttpClient
}
