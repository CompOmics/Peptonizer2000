use serde::Serialize;

/// A trait for performing HTTP POST requests with serialized payloads.
pub trait HttpClient {

    /// Performs an HTTP POST request with a JSON payload.
    ///
    /// # Arguments
    /// * `url` - Target URL for the request.
    /// * `batch` - Payload to serialize and send as JSON.
    ///
    /// # Returns
    /// Response body as a string on success.
    ///
    /// # Errors
    /// Returns an error string if the request fails or the response cannot be processed.
    fn perform_post_request<T: Serialize>(&self, url: String, batch: &T) -> Result<String, Box<dyn std::error::Error>>;
}

#[cfg(target_arch = "wasm32")]
pub struct WasmHttpClient;

#[cfg(not(target_arch = "wasm32"))]
pub struct PyHttpClient;

#[cfg(target_arch = "wasm32")]
impl HttpClient for WasmHttpClient {

    /// Sends an HTTP POST request in a WebAssembly environment using `XmlHttpRequest`.
    ///
    /// # Arguments
    /// * `url` - Target URL for the request.
    /// * `payload` - Payload to serialize and send as JSON.
    ///
    /// # Returns
    /// Response body as a string on success.
    ///
    /// # Errors
    /// Returns an error string if request setup, transmission, or response parsing fails.
    fn perform_post_request<T: Serialize>(&self, url: String, payload: &T) -> Result<String, Box<dyn std::error::Error>> {
        use web_sys::{XmlHttpRequest};

        let payload_json = serde_json::to_string(payload)?;

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

        Err(format!("Status code {}", status).into())
    }
}

#[cfg(not(target_arch = "wasm32"))]
impl HttpClient for PyHttpClient {

    /// Sends an HTTP POST request in a native environment using `reqwest` and Tokio.
    ///
    /// # Arguments
    /// * `url` - Target URL for the request.
    /// * `payload` - Payload to serialize and send as JSON.
    ///
    /// # Returns
    /// Response body as a string on success.
    ///
    /// # Errors
    /// Returns an error string if the request fails or the response cannot be processed.
    fn perform_post_request<T: Serialize>(&self, url: String, payload: &T) -> Result<String, Box<dyn std::error::Error>> {
        use reqwest::Client;
        use tokio::runtime::Runtime;

        // Create a Tokio runtime for async execution
        let rt = Runtime::new()?;

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
            Err(e) => Err(format!("HTTP POST request failed: {}", e).into()),
        }
    }
}

#[cfg(target_arch = "wasm32")]
/// Creates a new `HttpClient` for WebAssembly targets.
pub fn create_http_client() -> impl HttpClient {
    WasmHttpClient
}

#[cfg(not(target_arch = "wasm32"))]
/// Creates a new `HttpClient` for native targets.
pub fn create_http_client() -> impl HttpClient {
    PyHttpClient
}


#[cfg(test)]
mod tests {
    use super::*;
    use serde::Serialize;
    use serde_json::json;

    #[derive(Serialize)]
    struct DummyPayload {
        message: String,
    }

    #[test]
    fn test_invalid_url_returns_error() {
        let client = PyHttpClient;
        let payload = DummyPayload { message: "hello".to_string() };

        let result = client.perform_post_request("http://invalid_url".to_string(), &payload);
        assert!(result.is_err());
    }

    #[test]
    fn test_real_http_post() {
        let client = PyHttpClient;
        let payload = json!({ "foo": "bar" });

        let result = client.perform_post_request("https://api.unipept.ugent.be".to_string(), &payload);
        assert!(result.is_ok());
    }
}