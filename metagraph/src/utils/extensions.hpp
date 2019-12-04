#ifndef __EXTENSIONS_HPP__
#define __EXTENSIONS_HPP__

#include <functional>
#include <memory>
#include <assert.h>

#include "utils/string_utils.hpp"


namespace utils {

template <class T>
class Extendable {
  protected:
    class Extension {
      public:
        virtual ~Extension() {}
        virtual bool load(const std::string &filename_base) = 0;
        virtual void serialize(const std::string &filename_base) const = 0;
        virtual bool is_compatible(const T &obj, bool verbose = true) const = 0;
    };

  public:

    // TODO: improve interface: either prohibit or support
    //       properly multiple extensions of the same type.
    //       Use Extension::file_extension() or enum.
    void add_extension(std::shared_ptr<Extension> extension) {
        assert(extension.get());
        extensions_.push_back(extension);
    }

    template <class ExtensionSubtype>
    std::shared_ptr<ExtensionSubtype> get_extension() const {
        static_assert(std::is_base_of<Extension, ExtensionSubtype>::value);
        for (auto extension : extensions_) {
            if (auto match = std::dynamic_pointer_cast<ExtensionSubtype>(extension))
                return match;
        }
        return nullptr;
    }

    template <class ExtensionSubtype>
    void remove_extension() {
        static_assert(std::is_base_of<Extension, ExtensionSubtype>::value);
        for (auto it = extensions_.begin(); it != extensions_.end(); ++it) {
            if (auto match = std::dynamic_pointer_cast<ExtensionSubtype>(*it)) {
                extensions_.erase(it);
                return;
            }
        }
    }

    template <class ExtensionSubtype>
    void for_each(std::function<void(ExtensionSubtype &extension)> callback) {
        static_assert(std::is_base_of<Extension, ExtensionSubtype>::value);
        for (auto extension : extensions_) {
            if (auto match = std::dynamic_pointer_cast<ExtensionSubtype>(extension))
                callback(*match);
        }
    };

    template <class ExtensionSubtype>
    std::shared_ptr<ExtensionSubtype> load_extension(const std::string &filename) {
        static_assert(std::is_base_of<Extension, ExtensionSubtype>::value);
        remove_extension<ExtensionSubtype>();
        auto extension = std::make_shared<ExtensionSubtype>();

        auto filename_base = utils::remove_suffix(filename, file_extension());

        if (!extension->load(filename_base + file_extension()))
            return nullptr;

        add_extension(extension);
        return extension;
    }

    void serialize_extensions(const std::string &filename) const {
        for (auto extension : extensions_) {
            extension->serialize(utils::remove_suffix(filename, file_extension()) + file_extension());
        }
    }

    virtual std::string file_extension() const = 0;

  protected:
    std::vector<std::shared_ptr<Extension>> extensions_;
};

} // namespace utils

#endif // __EXTENSIONS_HPP__
