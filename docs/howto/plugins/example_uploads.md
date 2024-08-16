# How to write an example upload

Example uploads can be used to add representative collections of data for your plugin. Example uploads are available for end-users in the *Uploads*-page under the *Add example uploads*-button. There users can instantiate an example upload with a click. This can be very useful for educational or demonstration purposes but also for testing.

This documentation shows you how to write a plugin entry point for an example upload. You should read the [documentation on getting started with plugins](./plugins.md) to have a basic understanding of how plugins and plugin entry points work in the NOMAD ecosystem.

## Getting started

You can use our [template repository](https://github.com/FAIRmat-NFDI/nomad-plugin-template) to create an initial structure for a plugin containing an example upload. The relevant part of the repository layout will look something like this:

```txt
nomad-example
   ├── src
   │   ├── nomad_example
   │   │   ├── example_uploads
   │   │   │   ├── getting_started
   │   │   │   ├── __init__.py
   ├── LICENSE.txt
   ├── README.md
   ├── MANIFEST.in
   └── pyproject.toml
```

See the documentation on [plugin development guidelines](./plugins.md#plugin-development-guidelines) for more details on the best development practices for plugins, including linting, testing and documenting.

## Example upload entry point

The entry point is an instance of a `ExampleUploadEntryPoint` or its subclass. It defines basic information about your example upload and is used to automatically load the associated data into a NOMAD distribution. The entry point should be defined in `*/example_uploads/__init__.py` like this:

```python
from nomad.config.models.plugins import ExampleUploadEntryPoint

myexampleupload = ExampleUploadEntryPoint(
    title = 'My Example Upload',
    category = 'Examples',
    description = 'Description of this example upload.',
    path='example_uploads/getting_started
)
```

The default method for including the upload data is to place it in the plugin repository and use the `path` field to specify the location with respect to the package root. You can learn more about different data loading options in the next section. In the reference you can also see all of the available [configuration options for a `ExampleUploadEntryPoint`](../../reference/plugins.md#exampleuploadentrypoint).

The entry point instance should then be added to the `[project.entry-points.'nomad.plugin']` table in `pyproject.toml` in order for the example upload to be automatically detected:

```toml
[project.entry-points.'nomad.plugin']
myexampleupload = "nomad_example.example_uploads:myexampleupload"
```

## Including data in an example upload

There are three main ways to include data in an example upload:

1. Data stored directly in the plugin package using `path`:

    This is the default method that assumes you simply store the data under a path in the plugin source code. This is very convenient if you have relative small example data and wish to track this in version control. The path should be given relative to the package installation location (`src/<package-name>`), and you should ensure that the data is distributed with your Python package. Distribution of additional data files in Python packages is controlled with the `MANIFEST.in` file. If you create a plugin with our [template](https://github.com/FAIRmat-NFDI/nomad-plugin-template), the `src/<package-name>/example_uploads` folder is included automatically in `MANIFEST.in`. If you later add an example upload entry point to your plugin, remember to include the folder by adding the following line to `MANIFEST.in`:

    ```sh
    graft src/<package-name>/<path>
    ```

2. Data retrieved online during app startup using `url`:

    If your example uploads are very large (>100MB), storing them in Git may become unpractical. In order to deal with larger uploads, they can be stored in a separate online service. To load such external resources, you can specify a `url` parameter to activate online data retrieval. This will retrieve the large online file once upon the first app launch and then cache it for later use:

    ```python
    from nomad.config.models.plugins import ExampleUloadEntryPoint

    myexampleupload = ExampleUploadEntryPoint(
        name = 'MyExampleUpload',
        description = 'My custom example upload.',
        url='http://my_large_file_address.zip
    )
    ```

    Note that if the online file changes, you will need to remove the cached file for the new version to be retrieved. You can find the cached file in the package installation location, under folder `example_uploads`.

3. Data retrieved with a custom method:

    If the above options do not suite your use case, you can also override the `load`-method of `ExampleUploadEntryPoint` to perform completely custom data loading logic. Note that the loaded data should be saved in the package installation directory in order to be accessible. Check the default `load` function for more details.

    ```python
    from pydantic import Field
    from nomad.config.models.plugins import ExampleUploadEntryPoint


    class MyExampleUploadEntryPoint(ExampleUploadEntryPoint):

        def load(self):
            """Add your custom loading logic here."""
            ...
    ```