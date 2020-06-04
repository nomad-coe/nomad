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
    "rules": {
        "space-before-function-paren": ["error", "never"],
        "camelcase": [0]
    },
    "settings": {
        "react": {
            "version": "detect"
        }
    }
}