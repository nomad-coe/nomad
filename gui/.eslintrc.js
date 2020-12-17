module.exports = {
    "extends": [
        "standard",
        "plugin:react/recommended"
    ],
    "parser": "babel-eslint",
    "parserOptions": {
        "ecmaFeatures": {
            "jsx": true,
            "modules": true
        }
    },
    "globals": {
        "fetch": false,
        "browser": true
    },
    "plugins": [
        "react", "react-hooks", "testing-library", "jest"
    ],
    "rules": {
        "space-before-function-paren": ["error", {
            "asyncArrow": "always",
            "named": "never",
            "anonymous": "never"
        }],
        "camelcase": [0],
        "react-hooks/rules-of-hooks": "error",
        "react-hooks/exhaustive-deps": "warn",
        "react/display-name": [0]
    },
    "settings": {
        "react": {
            "version": "detect"
        }
    },
    "overrides": [
        {
            "files": [
                "**/*.spec.js",
                "**/*.spec.jsx"
            ],
            "env": {
                "jest": true
            }
        }
    ]
}